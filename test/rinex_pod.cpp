#include "beacon_tbl.hpp"
#include "datetime/datetime_write.hpp"
#include "datetime/utcdates.hpp"
#include "doris_observation_equations.hpp"
#include "doris_rinex.hpp"
#include "doris_utils.hpp"
#include "filters/filters.hpp"
#include "geodesy/geoconst.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iers2010.hpp"
#include "iers2010/iersc.hpp"
#include "iers2010/tropo.hpp"
#include "integrators.hpp"
#include "satellites/jason3.hpp"
#include "satellites/jason3_quaternions.hpp"
#include "sp3/sp3.hpp"
#include "sp3/sv_interpolate.hpp" /* debug mode */
#include "var_utils.hpp"
#include <cstdio>
#include <cassert>
#include "astrodynamics.hpp"

constexpr const int Np = 1;

// constexpr const double EleCutOff = 10e0; // elevation cut-off angle, [deg]

// max time difference between two observations to mark a new arc pass [sec]
constexpr const long max_sec_for_new_arc = 5 * 60L;

// default, i.e. a-priori value for Dfe / feN (relative frequency offset for
// emitter)
constexpr const double DfefeN_apriori = 0e0;
// respective default std. deviation
constexpr const double DfefeN_apriori_stddev = 1e0;

// a priori value for Drag Coefficient and respective sigma
constexpr const double DragCoef = 2e0;
constexpr const double DragCoefSigma = .05e0;

// Standard gravitational parameters for Sun and Moon in [km^3 / sec^2]
double GMSun, GMMoon;

// usually using these datetimes ...
using Datetime = dso::datetime<dso::nanoseconds>;

double crange(const Eigen::Matrix<double, 3, 1> &rbeacon_ecef,
              const Eigen::Matrix<double, 3, 1> &rsat_ecef,
              [[maybe_unused]] const dso::datetime<dso::nanoseconds> &ttai,
              Eigen::Matrix<double, 3, 1> &rbeacon_ecef_corr) noexcept {
  // return (rsat_ecef - rbeacon_ecef).norm()
  //   + iers2010::OmegaEarth * (rbeacon_ecef.dot(rsat_eci))/iers2010::C;
  Eigen::Matrix<double, 3, 1> newpos = rbeacon_ecef;
  rbeacon_ecef_corr = Eigen::Matrix<double, 3, 1>::Zero();
  int iteration = 0;
  constexpr const int max_iterations = 50;

  while ((rbeacon_ecef_corr - newpos).cwiseAbs().maxCoeff() > 1e-6 &&
         ++iteration < max_iterations) {
    rbeacon_ecef_corr = newpos;
    // t_emmision = t_reception - dt
    double r_approx = (rsat_ecef - rbeacon_ecef_corr).norm(); // range in [m]
    double dt = r_approx / iers2010::C;                       // [sec]
    newpos = Eigen::AngleAxisd(-iers2010::OmegaEarth * dt,
                               -Eigen::Vector3d::UnitZ()) *
             rbeacon_ecef;

    //Eigen::Matrix<double, 3, 1> diff = rbeacon_ecef_corr - newpos;
    //printf("\titeration, dt=%+.9f [sec] and DX=[%.9f, %.9f, %.9f]\n", dt,
    //       diff(0), diff(1), diff(2));
    //diff = newpos - rbeacon_ecef;
    //printf("\t                                  DX=[%.9f, %.9f, %.9f]\n",
    //       diff(0), diff(1), diff(2));
  }

  assert(iteration < max_iterations);

  return (rsat_ecef - newpos).norm();
}

// @see https://www.johndcook.com/blog/standard_deviation/
struct RunningStatistics {
  double _mean{0e0}, _var{0e0};
  long unsigned k{1};

  void update(double x) noexcept {
    double new_mean = _mean + (x - _mean) / (double)k;
    _var += (x - _mean) * (x - new_mean);
    _mean = new_mean;
    ++k;
  }

  long unsigned count() const noexcept {return k-1;}
  double mean() const noexcept { return _mean; }
  double variance() const noexcept { return _var / (k - 1); }
  double stddev() const noexcept { return std::sqrt(variance()); }
};

struct TropoDetails {
  double Lhz, mfh;
  double Lwz, mfw;
  double sum() const noexcept { return Lhz * mfh + Lwz * mfw; }
  double sum(double Lwz_) const noexcept { return Lhz * mfh + Lwz_ * mfw; }
  void dump(double Lwz_) const noexcept {
    printf("##TRP Tropo details: %.6f * %.6f + (%.6f or %.6f) * %.6f\n", Lhz, mfh, Lwz, Lwz_, mfw);
  }
};

struct SatBeacon {
  char id3c[3];                  ///< beacons internal (3-char) id
  int count;                     ///< internally used integer id
  Datetime ttai, tproper;        ///< time
  double Ls1, Lu2;               ///< cycles on L1 (S1)
  double Diono;                  ///< ionospheric corrections [cycles]
  double Drel;                   ///< relativity correction
  TropoDetails Dtropo;           ///< tropospheric correction
  Eigen::Matrix<double, 3, 1> s; ///< beacon-satellite vector in topocentric rf
  double crho;                   ///< beacon-satellite range, corrected for Earth rotation (Sagnac & aberration)
  int arcnr{0};
  int reinitialize{0}; ///< measurement marked as discontinuity

  SatBeacon(const char *id_, int count_, const Datetime &ttai_,
            const Datetime &tproper_, double L1, double L2, double Diono_,
            const TropoDetails &Dtropo_, double Drel_,
            const Eigen::Matrix<double, 3, 1> &s_, double _crho) noexcept
      : count(count_), ttai(ttai_), tproper(tproper_), Ls1(L1), Lu2(L2),
        Diono(Diono_), Drel(Drel_), Dtropo(Dtropo_), s(s_), crho(_crho) {
    // std::strncpy(id3c, id_, 3); // fuck the warning!
    std::memcpy(id3c, id_, sizeof(char) * 3);
  }

  SatBeacon &operator=(const SatBeacon &sb) noexcept {
    std::strncpy(id3c, sb.id3c, 3);
    count = sb.count;
    ttai = sb.ttai;
    tproper = sb.tproper;
    Ls1 = sb.Ls1;
    Lu2 = sb.Lu2;
    Diono = sb.Diono;
    Dtropo = sb.Dtropo;
    Drel = sb.Drel;
    s = sb.s;
    crho = sb.crho;
    return *this;
  }

  void update(const Datetime &ttai_, const Datetime &tproper_, double L1_,
              double L2_, double Diono_, const TropoDetails &Dtropo_,
              double Drel_, const Eigen::Matrix<double, 3, 1> &s_, double _crho) noexcept {
    ttai = ttai_;
    tproper = tproper_;
    Ls1 = L1_;
    Lu2 = L2_;
    Diono = Diono_;
    Dtropo = Dtropo_;
    Drel = Drel_;
    s = s_;
    crho = _crho;
  }

  double rho() const noexcept { return /*s.norm();*/crho; }
};

// hold satellite state & time
struct SatelliteState {
  // Current datetime in TAI
  double mjd_tai;
  // state vector at t=tai in ECEF, at DORIS receiver RP (Iono-Free)
  Eigen::Matrix<double, 6, 1> state;
  // state transition matrix at t=tai
  Eigen::Matrix<double, 6, 6> Phi;
  Eigen::Matrix<double, 6, Np> S;
  // Satellite's center of gravity w.r.t the satellite-fixed frame
  Eigen::Matrix<double, 3, 1> *cog_sf{nullptr};
  // Satellite's 2GHz receiver ARP w.r.t the satellite-fixed frame
  Eigen::Matrix<double, 3, 1> *arp_sf{nullptr};
  // Body-quaternion hunter
  dso::JasonQuaternionHunter *qhunt{nullptr};

  // transform state vector state from ECEF to GCRF for mjd_tai
  Eigen::Matrix<double, 6, 1>
  celestial(const dso::EopLookUpTable &elut) const noexcept {
    Eigen::Matrix<double, 3, 3> rc2i, rpom;
    double era, lod;
    assert(!gcrs2itrs(mjd_tai, elut, rc2i, era, rpom, lod));
    return dso::yter2cel(state, rc2i, era, lod, rpom);
  }

  // Vector to go from receiver ARP to satellite's CoG, in GCRS
  // aka dx = CoG - ARP in GCRF
  Eigen::Matrix<double, 3, 1>
  eccentricity([[maybe_unused]] const Eigen::Quaternion<double> &q) noexcept {
    // assuming that the quaternion acts as:
    // X_sat-fixed = Q * X_sat-gcrs
    return (cog_sf) ? (q.conjugate().normalized() * (*cog_sf - *arp_sf))
                    : (Eigen::Matrix<double, 3, 1>::Zero());
  }

  int integrate(double mjd_target, dso::SGOde &integrator) noexcept {
    // count calls; this is only needed because the first call should
    // consider the satellite coordinates as CoM coordinates, and not
    // apply eccentricity (ARP to  CoM)
    static int call_nr = 0;

    // Vector containing state + variational equations size: 6 + 6x(6+Np)
    // Stored in plain format, one column at a time
    // Ref. Frame: inertial
    Eigen::Matrix<double, 6 + 6 * 6 + 6 * Np, 1> yPhi =
        Eigen::Matrix<double, 6 + 6 * 6 + 6 * Np, 1>::Zero();

    // transform state from ECEF to GCRF (satellite ARP)
    yPhi.block<6, 1>(0, 0) = celestial(integrator.params->eopLUT);

#ifndef NO_ATTITUDE
    // if needed, go from antenna RP to CoG (this is not needed in the first
    // call, since we integrate starting with the sp3 state which is w.r.t.
    // the satellite's CoG)
    if (call_nr) {
      Eigen::Quaternion<double> q;
      if (qhunt->get_at(mjd_tai, q)) {
        fprintf(stderr, "ERROR Failed to find quaternion for datetime\n");
        assert(false);
      }
      // GCRF: from satellite ARP to CoG
      yPhi.block<3, 1>(0, 0) += eccentricity(q);
    }
#endif

    // Initial condition for state transition matrix Φ(t0,t0) = I
    // set the state transition matrix to identity (initial condition)
    // Note that the matrix is stored in plain format, one column at a time
    // First column contains the state, skip it.
    {
      int k = 6;
      for (int col = 1; col < 6; col++)
        for (int row = 0; row < 6; row++)
          yPhi(k++) = (col - 1 == row) ? 1e0 : 0e0;
    }

    // Initial condition for state transition matrix S(t0) = I
    // last 6 * Np elements should be zero (already done in initialization of
    // the matrix)

    // t0 for variational equations (TAI)
    integrator.params->mjd_tai = mjd_tai;

    // target t for variational equations; seconds after t0
    double tout = (mjd_target - mjd_tai) * 86400e0;

    // set initial intergation flag
    integrator.flag() = 1;

    // keep solution here (celestial RF at tout)
    // As yPhi, this is in plain format, one column at a time
    Eigen::VectorXd sol(6 + 6 * 6 + 6 * Np);

    // integrate (in inertial RF), from 0+mjd_tai to tout+mjd_tai [sec]
    double tsec = 0e0;
    integrator.de(tsec, tout, yPhi, sol);

    // output epoch as datetime
    double tout_mjd = integrator.params->mjd_tai + tsec / 86400e0;

    // let's see were we are at
    if (std::abs(tout_mjd - mjd_target) > 1e-12) {
      fprintf(stderr,
              "ERROR wanted integration to %.9f and got up to %.9f, that is "
              "%.9f sec apart!\n",
              mjd_target, tout_mjd, (mjd_target - mjd_tai) * 86400e0);
      return 1;
    }

    // everything seems ok, update state and time
    mjd_tai = tout_mjd;

#ifndef NO_ATTITUDE
    // if needed, go from CoM to antenna ARP (sol must be in GCRS)
    {
      Eigen::Quaternion<double> q;
      if (qhunt->get_at(mjd_tai, q)) {
        fprintf(stderr, "ERROR Failed to find quaternion for datetime\n");
        assert(false);
      }
      sol.block<3, 1>(0, 0) -= eccentricity(q);
    }
#endif

    // Transform state back to ECEF (from GCRF), wrt the satellite ARP
    Eigen::Matrix<double, 3, 3> rc2i, rpom;
    double era, lod;
    assert(
        !gcrs2itrs(mjd_tai, integrator.params->eopLUT, rc2i, era, rpom, lod));
    state = dso::ycel2ter(sol.block<6, 1>(0, 0), rc2i, era, lod, rpom);

    // assign Phi matrix (6x6)
    for (int i = 0; i < 6; i++) {
      Phi.col(i) = yPhi.block<6, 1>(6 * (i + 1), 0);
    }

    // assign S matrix (6xNp)
    if (Np) {
      for (int i=0; i<Np; i++) {
        S.col(i) = yPhi.block<6, 1>(6 * (i + 1) + 6*6, 0);
      }
    }

    ++call_nr;
    return 0;
  }
};

int relativity_corrections(const Eigen::Matrix<double, 6, 1> &sv_state,
                           const Eigen::Matrix<double, 3, 1> &rbeacon,
                           double Re, double GM, double J2, double &Drel_c,
                           double &Drel_r) noexcept;

int prepare_beacon_coordinates(
    std::vector<dso::BeaconCoordinates> &beaconCrdVec,
    const char *beacon_info_tbl,
    const dso::datetime<dso::nanoseconds> &t) noexcept;

Eigen::Matrix<double, 3, 1>
beacon_arp2ion(const Eigen::Matrix<double, 3, 1> &bxyz_arp,
               const dso::BeaconStation &beacon) noexcept;

// get the l1, l2 and f indexes off from a RINEX file (instance)
int get_rinex_indexes(const dso::DorisObsRinex &rnx, int &l1, int &l2, int &f,
                      int &w1_idx, int &w2_idx) noexcept;

Eigen::Matrix<double, 3, 1>
beacon_coordinates(const char *_4charid,
                   const std::vector<dso::BeaconCoordinates> &crdVec) noexcept;

int get_tropo(const dso::datetime<dso::nanoseconds> &t,
              const Eigen::Matrix<double, 3, 1> &bxyz, double zd,
              const dso::Gpt3Grid &grid, TropoDetails &Dtrop) noexcept;
int get_tropo_vmf(const char *site, const dso::datetime<dso::nanoseconds> &t,
              const Eigen::Matrix<double, 3, 1> &bxyz, double zd,
              dso::SiteVMF3Feed& feed, TropoDetails &Dtrop) noexcept;

Eigen::Matrix<double, 3, 1>
range_rate(const Eigen::Matrix<double, 3, 1> &s_t0,
           const Eigen::Matrix<double, 3, 1> &rsta,
           const Eigen::Matrix<double, 3, 1> &rsat, double dt,
           double &rr_computed, Eigen::Matrix<double, 6, 1> &drrdrv) noexcept;

std::vector<const char *>
beaconcrs2cchar_vec(const std::vector<dso::BeaconCoordinates> &bcv) noexcept {
  std::vector<const char *> site_names;
  site_names.reserve(bcv.size());
  for (auto it = bcv.begin(); it != bcv.end(); ++it) {
    site_names.push_back(it->id);
  }
  return site_names;
}

int main(int argc, char *argv[]) {
  // check input
  if (argc != 2) {
    fprintf(stderr, "USAGE: %s [YAML CONFIG]\n", argv[0]);
    return 1;
  }

  // resolve the yaml config file and get the root node
  const YAML::Node config = YAML::LoadFile(argv[1]);
  char buf[256];
  int error;

  // Elevation Cut-Off Angle in [degrees]
  // -------------------------------------------------------------------------
  double EleCutOff;
  if (dso::get_yaml_value_depth2<double>(config, "filtering", "elevation-cut-off", EleCutOff)) {
    fprintf(stderr, "ERROR failed to find spk kernel\n");
    return 1;
  }
  printf("# Using elevation cut-off angle = %.3f\n", EleCutOff);

  // Load CSPICE/NAIF Kernels (2/3)
  // -------------------------------------------------------------------------
  if (dso::get_yaml_value_depth2(config, "naif-kernels", "spk", buf)) {
    fprintf(stderr, "ERROR failed to find spk kernel\n");
    return 1;
  }
  dso::cspice::load_if_unloaded_spk(buf);
  if (dso::get_yaml_value_depth2(config, "naif-kernels", "lsk", buf)) {
    fprintf(stderr, "ERROR failed to find spk kernel\n");
    return 1;
  }
  dso::cspice::load_if_unloaded_lsk(buf);

  // DORIS RINEX instance
  // -------------------------------------------------------------------------
  if (dso::get_yaml_value_depth2(config, "data", "doris-rinex", buf)) {
    fprintf(stderr, "ERROR. Failed parsing data/rinex file from YAML %s\n",
            argv[1]);
    return 1;
  }

  // Fit a linear Model for the Relative Receiver Offset Values (use proper 
  // time for this).
  dso::PolynomialModel<dso::datetime<dso::nanoseconds>> rfo_fit(1);
  {
    if (dso::fit_relative_frequency_offset(buf, rfo_fit, false, 4)) {
      fprintf(stderr, "ERROR! Failed to fit RFO values!\n");
      return 2;
    }
  }

  // Construct the Doris RINEX instance rnx
  dso::DorisObsRinex rnx(buf);

  // Initial Orbit
  // -------------------------------------------------------------------------
  // Intial satellite state, get it from the Sp3 using the RINEX's time of
  // first obs
  SatelliteState svState;
  dso::get_yaml_value_depth2(config, "data", "sp3", buf);
  dso::Sp3c sp3(buf);
  dso::Sp3Iterator sp3_iterator(sp3);

  // EOP Look Up Table
  // Parse the input EOP data file to create an EopLookUpTable eop_lut
  // -------------------------------------------------------------------------
  dso::EopLookUpTable eop_lut;
  if (dso::get_yaml_value_depth2(config, "eop-info", "eop-file", buf)) {
    fprintf(stderr,
            "ERROR. Failed parsing eop-info/eop-file file from YAML %s\n",
            argv[1]);
    return 1;
  } else {
    const int ref_mjd = rnx.ref_datetime().as_mjd();
    const int start = ref_mjd - 5;
    const int end = ref_mjd + 6;
    // parse C04 EOPs and convert time-stamps to TT (not UTC)
    if (parse_iers_C04(buf, start, end, eop_lut)) {
      fprintf(stderr, "ERROR. Failed collecting EOP data\n");
      return 1;
    }
    // regularize ERP (DUT and LOD)
    eop_lut.regularize();
  }

  // Gravity
  // -------------------------------------------------------------------------
  // parse degree, order and the requested gravity model into a
  // HarmonicCoeffs instance. Note that to compute potential we will need
  // Lagrange polynomials (later on)
  int degree, order;
  error = 0;
  error = dso::get_yaml_value_depth2<int>(config, "gravity", "degree", degree);
  error += dso::get_yaml_value_depth2<int>(config, "gravity", "order", order);
  dso::HarmonicCoeffs harmonics(degree);
  if (!error)
    error = dso::get_yaml_value_depth2(config, "gravity", "model", buf);
  if (!error)
    error = dso::parse_gravity_model(buf, degree, order, rnx.ref_datetime(),
                                     harmonics, true);
  if (error) {
    fprintf(stderr, "ERROR Failed handling gravity field model!\n");
    return 1;
  }

  // Station/Beacon coordinates
  // -------------------------------------------------------------------------
  // Get beacon coordinates from sinex file and extrapolate to RINEX ref. time
  // Result coordinates per beacon are stored in the beaconCrdVec vector.
  // Note that these poition vectors are w.r.t the beacon/antenna reference
  // point. When in actual processing, this has to be changed, if we are
  // considering iono-free analysis
  // -------------------------------------------------------------------------
  if (dso::get_yaml_value_depth2(config, "reference-frame",
                                 "station-coordinates", buf)) {
    fprintf(stderr,
            "ERROR. Failed parsing reference-frame/station-coordinates file "
            "from YAML %s\n",
            argv[1]);
    return 1;
  }
  std::vector<dso::BeaconCoordinates> beaconCrdVec;
  beaconCrdVec.reserve(rnx.stations().size());
  // coordinates of beacons (on ground), ECEF/pdop
  if (extrapolate_sinex_coordinates(buf, rnx.stations(), rnx.ref_datetime(),
                                    beaconCrdVec, false)) {
    fprintf(stderr,
            "ERROR. Failed extracting/extrapolating beacon coordinates\n");
    return 1;
  }
  // coordinates of beacons (on antenna RP), ECEF/pdop
  dso::get_yaml_value_depth2(config, "data", "beacon-information", buf);
  if (prepare_beacon_coordinates(beaconCrdVec, buf, rnx.ref_datetime())) {
    fprintf(stderr, "ERROR Failed applying eccentricities to beacons!\n");
    return 1;
  }

  // Troposphere
  // -------------------------------------------------------------------------
  // read-in the grid file. That is all for now
  /* for gpt3/vmf3 
  dso::get_yaml_value_depth3(config, "troposphere", "gpt3", "grid", buf);
  dso::Gpt3Grid gpt3_grid(buf);
  */
  dso::get_yaml_value_depth3(config, "troposphere", "vmf3", "grid", buf);
  auto site_names_vec = beaconcrs2cchar_vec(beaconCrdVec);
  dso::SiteVMF3Feed feed(buf, site_names_vec);

  // Attitude Information
  dso::get_yaml_value_depth2(config, "attitude", "body-quaternion", buf);
  dso::JasonQuaternionHunter qhunt(buf);

  // Previous Observation relating Beacon/Satellite (to compute Ndop)
  // -------------------------------------------------------------------------
  std::vector<SatBeacon> prevec;
  prevec.reserve(beaconCrdVec.size());

  // On-board receiver eccentricity, in the satellite-fixed frame
  // -------------------------------------------------------------------------
  if (dso::get_yaml_value_depth2(config, "attitude", "mass-cog", buf)) {
    fprintf(stderr, "ERROR Failed locating Mass and CoG information file\n");
    return 1;
  }
  Eigen::Matrix<double, 3, 1> sat_cog, l3_pco;
  double sat_mass;
  {
    // get satellite CoG coordinates in the satellite-fixed RF (with
    // corrections)
    assert(!dso::SatelliteInfo<dso::SATELLITE::Jason3>::mass_cog(
        rnx.ref_datetime(), buf, sat_mass, sat_cog));
    // get satellite ARP coordinates in the satellite-fixed RF
    Eigen::Matrix<double, 3, 1> l1_pco, l2_pco;
    dso::SatelliteInfo<dso::SATELLITE::Jason3>::pco(l1_pco, l2_pco);
    // compute the iono-free phase center
    l3_pco = l1_pco + (l2_pco - l1_pco) / (dso::GAMMA_FACTOR - 1e0);
  }
  svState.cog_sf = &sat_cog;
  svState.arp_sf = &l3_pco;
  svState.qhunt = &qhunt;

  // Setup Integration Parameters for Orbit Integration
  // -------------------------------------------------------------------------
  // We will need the pck (SPICE) kernel for gravitational parameters of Sun
  // and Moon
  if (dso::get_yaml_value_depth2(config, "naif-kernels", "pck", buf)) {
    fprintf(stderr, "ERROR Failed locating NAIF pck kernel\n");
    return 1;
  }
  dso::IntegrationParameters IntegrationParams(degree, order, eop_lut,
                                               harmonics, buf);
  IntegrationParams.macromodel =
      dso::MacroModel<dso::SATELLITE::Jason3>::mmcomponents;
  IntegrationParams.numMacroModelComponents =
      dso::MacroModel<dso::SATELLITE::Jason3>::NumPlates;
  IntegrationParams.qhunt = &qhunt;
  IntegrationParams.SatMass = &sat_mass;

  dso::Nrlmsise00 nrlmsise00;
  IntegrationParams.nrlmsise00 = &nrlmsise00;

  dso::get_yaml_value_depth3(config, "force-model", "atmospheric-drag",
                             "atmo-data-csv", buf);
  dso::modified_julian_day utc_mjd_fo;
  const double utc_fday = dso::tai2utc(rnx.time_of_first_obs(), utc_mjd_fo);
  dso::nrlmsise00::InParams<
      dso::nrlmsise00::detail::FluxDataFeedType::ST_CSV_SW>
      atm_data_feed(buf, utc_mjd_fo, utc_fday * 86400e0);
  atm_data_feed.params_.set_switches_on();
  atm_data_feed.params_.use_aparray();
  atm_data_feed.params_.meters_on();
  IntegrationParams.AtmDataFeed = &atm_data_feed;

  // Orbit Integrator
  // -------------------------------------------------------------------------
  // Setup an integrator, to extrapolate orbit with:
  // 1. Relative accuracy 1e-12
  // 2. Absolute accuracy 1e-12
  // 3. Num of Equations: 6 for state and 6*6 for variational equations
  dso::SGOde Integrator(dso::VariationalEquations, 6 + 6 * 6 + 6 * Np, 1e-12,
                        1e-12, &IntegrationParams);

  // get the (RINEX) indexes for the observables we want
  int l1i, l2i, fi, w1i, w2i;
  if (get_rinex_indexes(rnx, l1i, l2i, fi, w1i, w2i))
    return 1;

  // Extended Kalman Filter instance
  // -------------------------------------------------------------------------
  //const int NumParams =
  //    6                        ///< satellite state
  //    + rnx.stations().size()  ///< beacon relative frequency offset
  //    + rnx.stations().size(); ///< wet tropo path dealy (zenith)
  // dso::ExtendedKalmanFilter<dso::nanoseconds> Filter(NumParams);
#ifdef RANDOM_RFO
  dso::EkfFilter<Np, dso::BeaconClockModel::None> Filter(rnx.stations().size(),
                                                    rnx.time_of_first_obs());
#else
  printf("Note: We are estimating linear relative frequency offsets!\n");
  dso::EkfFilter<Np, dso::BeaconClockModel::Linear> Filter(rnx.stations().size(),
                                                    rnx.time_of_first_obs());
#endif
  const int NumParams = Filter.num_params();

  // initialize the Kalman filter
#ifdef RANDOM_RFO
  Filter.set_rfoff_apriori(DfefeN_apriori);
  Filter.set_rfoff_apriori_sigma(DfefeN_apriori_stddev * DfefeN_apriori_stddev);
#else
  Filter.set_rfoff_apriori(DfefeN_apriori, 0e0);
  Filter.set_rfoff_apriori_sigma(DfefeN_apriori_stddev * DfefeN_apriori_stddev,
                                 DfefeN_apriori_stddev * DfefeN_apriori_stddev);
#endif
  Filter.set_tropo_apriori_sigma(.5e0 * .5e0);
  Filter.set_drag_coef_apriori(DragCoef);
  Filter.set_drag_coef_apriori_sigma(DragCoefSigma * DragCoefSigma);

  // Default observation sigma for a range-rate observable at zenith
  double obs_sigma;
  error = dso::get_yaml_value_depth2<double>(config, "filtering",
                                             "observation-sigma", obs_sigma);
  assert(obs_sigma > 0e0 && obs_sigma < 1e2 && (!error));
  printf("## Note: default observation sigma value: %.3f\n", obs_sigma);

  // Important !!
  // set integration parameter estimates to point to the Filter
  // IntegrationParams.estimates = &Filter.x;
  IntegrationParams.drag_coef = &(Filter._ekf._ekf.x(Filter.drag_coef_index()));

  // Start RINEX data-block iteration
  // -------------------------------------------------------------------------
  // get an iterator to the RINEXs data blocks
  dso::RinexDataBlockIterator it(&rnx);

  // Some variables ...
  const double J2 = harmonics.J2();
  const double GM = harmonics.GM();
  const double Re = harmonics.Re();

  // Running average and std. deviation of O-C values
  RunningStatistics rstats;

  error = 0;
  [[maybe_unused]] int dummy_counter = 0;
  // count beacons as we encounter them (above min. elevation)
  int receiver_count = 0;
  // count Ndop observations
  unsigned ndop_count = 0;
  unsigned ndop_count_rejected = 0;
  unsigned flaged_obs = 0; // number of observations with 'bad' flags
  unsigned num_obs = 0; // observation count, regardless if usable or not

  // report
  #ifndef NO_ATTITUDE
    printf("## Attitude information will be used!\n");
  #else
    printf("## No attitude information will be used!\n");
  #endif

  // for every new data block in the RINEX file (aka every epoch) ...
  while (!(error = it.next())) {

    // the current reference time for the L1 observation (corrected for
    // receiver clock offset). That is tl1 is approximately TAI.
    // aka proper-time to coordinate-time
    const auto tl1 = it.corrected_l1_epoch();

    // current proper time (aka tau)
    const auto tproper = it.proper_time();

    // a buffer to write datetime strings to
    char dtbuf[64];

    // integrate orbit to here (TAI)
    // first get reference state from sp3, for an epoch as close as possible
    if (!dummy_counter) {
      if (sp3_iterator.goto_epoch(tl1)) {
        fprintf(stderr, "ERROR Failed to get reference position from SP3\n");
        return 1;
      }
      svState.mjd_tai = sp3_iterator.current_time().as_mjd();
      svState.state(0) = sp3_iterator.data_block().state[0] * 1e3;
      svState.state(1) = sp3_iterator.data_block().state[1] * 1e3;
      svState.state(2) = sp3_iterator.data_block().state[2] * 1e3;
      svState.state(3) = sp3_iterator.data_block().state[4] * 1e-1;
      svState.state(4) = sp3_iterator.data_block().state[5] * 1e-1;
      svState.state(5) = sp3_iterator.data_block().state[6] * 1e-1;
      if (svState.integrate(tl1.as_mjd(), Integrator)) {
        fprintf(stderr, "ERROR. Failed to integrate orbit!\n");
        return 1;
      }
    } else {
      if (svState.integrate(tl1.as_mjd(), Integrator)) {
        fprintf(stderr, "ERROR. Failed to integrate orbit!\n");
        return 1;
      }
    }

    // We now have the state vector svState.state in ECEF coordinates wrt
    // the satellite antenna RP. The state vector referes to svState.tai_mjd.
    // We also have the state transition matrix ready
    
    // NEW::
    // GCRS-to-ITRS parameters for current TAI
    Eigen::Matrix<double, 3, 3> rc2i_now, rpom_now;
    double era_now, lod_now;
    assert(!gcrs2itrs(tl1.as_mjd(), eop_lut, rc2i_now, era_now, rpom_now,
                      lod_now));
    // return dso::yter2cel(state, rc2i, era, lod, rpom);

    // needed for tides ...
    // Various variables to be used, depending only on the current epoch
    //Eigen::Matrix<double, 3, 1> rmoon, rsun; // [km], ITRF
    //dso::datetime<dso::nanoseconds> tt(tl1);
    //double utc_fhours;
    //{
    //  tt.add_seconds(dso::nanoseconds(32184 * 1'000'000L));
    //  double vec[3];
    //  int terror = dso::moon_vector_cspice(tt, vec);
    //  rmoon = Eigen::Matrix<double, 3, 1>(vec); // [km] GCRF
    //  terror += dso::sun_vector_cspice(tt, vec);
    //  rsun = Eigen::Matrix<double, 3, 1>(vec); // [km] GCRF
    //  assert(!terror);
    //  Eigen::Matrix<double, 3, 3> rc2i, rpom;
    //  double era, lod;
    //  assert(!gcrs2itrs(tl1.as_mjd(), eop_lut, rc2i, era, rpom, lod));
    //  rmoon = dso::rcel2ter(rmoon, rc2i, era, rpom);
    //  rsun = dso::rcel2ter(rsun, rc2i, era, rpom);
    //  dso::modified_julian_day mjdi;
    //  utc_fhours = dso::tai2utc(tl1, mjdi) * 24e0;
    //}

    // Kalman filter time update:
    // estimates: the state part is updated using the result of integration;
    //            the rest of the (estimated) parameters remain the same
    // P matrix : Update using the state transition matrix (result of orbit
    //            integration)
    Eigen::MatrixXd PhiP = Eigen::MatrixXd::Identity(NumParams, NumParams);
    PhiP.block<6, 6>(0, 0) = svState.Phi;
    auto estimates = Filter.estimates();
    estimates.block<6, 1>(0, 0) = svState.state;
    Filter.time_update(tl1, estimates, PhiP);
    #ifdef FIX_ORBIT
    Filter.P().block<6,6>(0,0) = Eigen::Matrix<double,6,6>::Identity() * 1e-12;
    #endif

    ++dummy_counter;

    // iterate through the observation set (aka the various beacons with
    // observations for current epoch)
    auto beaconobs = it.cblock.begin();
    while (beaconobs != it.cblock.end()) {

      // increase observation count
      ++num_obs;

      // an iterator to the RINEX's stations vector (aka a BeaconStation)
      // matching the current beacon
      auto beacon_it = rnx.beacon_internal_id2BeaconStation(beaconobs->id());
      assert(beacon_it != rnx.stations().cend());

      // check flags
      if (beaconobs->m_values[l1i].m_flag1 ==
              '1' || // 2 GHz central frequency measurement
          beaconobs->m_values[l1i].m_flag2 ==
              '1' || // discontinuity of 2 GHZ measurement
          beaconobs->m_values[l2i].m_flag1 ==
              '1' || // 400 MHz central frequency measurement
          beaconobs->m_values[l2i].m_flag2 ==
              '1' || // discontinuity of 400 MHZ measurement
          beaconobs->m_values[w1i].m_flag1 == '1' || // station on restart mode
          beaconobs->m_values[w2i].m_flag1 == '1'    // station on restart mode
      ) {
        // Oops! flags not ok! we are skipping this observation. first,
        // check if the observation is in the middle of the arc. If yes,
        // then mark a discontinuity
        auto pprev_obs = std::find_if(
            prevec.begin(), prevec.end(), [&](const SatBeacon &sb) {
              return (sb.id3c[0] == beaconobs->id()[0] &&
                      sb.id3c[1] == beaconobs->id()[1] &&
                      sb.id3c[2] == beaconobs->id()[2]);
            });
        // mark as discontinuity and skip observation
        if (pprev_obs != prevec.end()) {
          pprev_obs->reinitialize = 1;
        }

        // increase observation count for flagged observations
        ++flaged_obs;

      } else {
        // flags ok, continue with processing ...

        // we are going to need the beacons ECEF coordinates (note that
        // these are antenna RP coordinates)
        Eigen::Matrix<double, 3, 1> bxyz_sta =
            beacon_coordinates(beacon_it->m_station_id, beaconCrdVec);
        
        // add tidal displacement (dehanttideinel)
        //if (1==2) {
        //  Eigen::Matrix<double, 3, 1> toff = iers2010::dehanttideinel_impl(
        //      bxyz_sta, rsun * 1e3, rmoon * 1e3, tt.jcenturies_sinceJ2000(),
        //      utc_fhours);
        //  bxyz_sta = bxyz_sta + toff;
        //}

        // Iono-Free phase center, ECEF (changeme)
        Eigen::Matrix<double, 3, 1> bxyz_ion =
            beacon_arp2ion(bxyz_sta, *beacon_it);

        // get azimouth [rad], elevation [rad] and geometric distance [m]
        // (beacon to satellite, Iono-Free PC at both ends)
        Eigen::Matrix<double, 3, 1> r_enu =
            dso::car2top<dso::ellipsoid::grs80>(
                bxyz_ion, svState.state.block<3, 1>(0, 0));
        double az, el;
        double rho = dso::top2dae(r_enu, az, el);

        // correct the satellite-beacon distance for Earth rotation
        {
          Eigen::Matrix<double, 3, 1> corrected_beacon;
          rho = crange(bxyz_ion, svState.state.block<3, 1>(0, 0), tl1,
                       corrected_beacon);
          bxyz_ion = corrected_beacon;
          // !! WARNING !!
          r_enu = bxyz_ion;
        }

        // only process observations to elevation > EleCutOff [deg]
        if (dso::rad2deg(el) < EleCutOff) {
          ;

        } else {
          // elevation ok, continue processing

          // is this the start of a new arc ?
          bool start_new_arc = false;

          // check if we already have a previous observation for this
          // beacon
          auto pprev_obs = std::find_if(
              prevec.begin(), prevec.end(), [&](const SatBeacon &sb) {
                return (sb.id3c[0] == beaconobs->id()[0] &&
                        sb.id3c[1] == beaconobs->id()[1] &&
                        sb.id3c[2] == beaconobs->id()[2]);
              });
          // if we do, check if we have to start a new arc. Check
          // performed in proper time
          if (pprev_obs != prevec.end()) {
            if (tproper.delta_sec(pprev_obs->tproper) >
                dso::nanoseconds(max_sec_for_new_arc *
                                 dso::nanoseconds::sec_factor<
                                     dso::nanoseconds::underlying_type>())) {
              start_new_arc = true;
            }
          }

          // nominal frequencies for the beacon: feN
          int k; // shift factor
          if (rnx.beacon_shift_factor(beaconobs->id(), k)) {
            fprintf(stderr, "Failed to find shift factor for beacon %.3s\n",
                    beaconobs->id());
            return 1;
          }
          double feN, fe2N; // [Hz]
          dso::beacon_nominal_frequency(k, feN, fe2N);

          // ionospheric correction (actually L_2GHz = L_2GHz + Diono)
          const double cDiono =
              dso::carrier_iono_correction(beaconobs->m_values[l1i].m_value,
                                           beaconobs->m_values[l2i].m_value);

          // tropospheric correction
          TropoDetails cDtropo;
          /*if (get_tropo(tl1, bxyz_ion, dso::DPI / 2e0 - el, gpt3_grid,
                        cDtropo)) {
            return 3;
          }*/
          char __siteid5[5] = {'\0'};
          std::memcpy(__siteid5, beacon_it->m_station_id, 4 * sizeof(char));
          if (!get_tropo_vmf(__siteid5, tl1, bxyz_ion, dso::DPI / 2e0 - el, feed,
                            cDtropo)) 
          {

          // relativistic correction
          double Drel_c, Drel_r;
          relativity_corrections(svState.state, bxyz_ion, Re, GM, J2, Drel_c,
                                 Drel_r);
          const double cDrel = Drel_c + Drel_r;

          // Are we going to compute Doppler count / observations
          // equation? We may not need to, if :
          // 1. this is the first observation of the beacon
          //    aka pprev_obs == prevec.end()
          // 2. We are starting a new arc (aka start_new_arc==1)
          // 3. Previous observation was marked
          if (pprev_obs == prevec.end() || pprev_obs->reinitialize ||
              start_new_arc) {
            if (pprev_obs != prevec.end()) {
              // not the first observation for the beacon; update last
              // observation info
              pprev_obs->update(tl1, tproper, beaconobs->m_values[l1i].m_value,
                                beaconobs->m_values[l2i].m_value, cDiono,
                                cDtropo, cDrel, r_enu, rho);
              // update arc number -- if needed
              if (start_new_arc) {
                pprev_obs->arcnr += 1;
                // update a-priori value for this beacon zenith wet tropo
                // delay [m]
                Filter.tropo_estimate(pprev_obs->count) = cDtropo.Lwz;
              }
              // passed discontinuity, observation ok
              if (pprev_obs->reinitialize)
                pprev_obs->reinitialize = 0;
            } else {
              // first observation for the beacon
              prevec.emplace_back(SatBeacon(beaconobs->id(), receiver_count,
                                            tl1, tproper,
                                            beaconobs->m_values[l1i].m_value,
                                            beaconobs->m_values[l2i].m_value,
                                            cDiono, cDtropo, cDrel, r_enu, rho));
              // update a-priori value for this beacon zenith wet tropo
              // delay [m]
              Filter.tropo_estimate(receiver_count) = cDtropo.Lwz;
              // for next beacon count
              ++receiver_count;
            }

          } else {
            // compute Ndop and observation equation
            
            // The number/count of the beacon in the Filter
            const int receiver_number = pprev_obs->count;

            // we need to find the true proper frequency of the receiver
            // (aka satellite), f_rT [Hz]
            const double frT = dso::DORIS_FREQ1_MHZ * 1e6 *
                               (1e0 + rfo_fit.value_at(tproper) * 1e-11);
            //const double frT = dso::DORIS_FREQ1_MHZ * 1e6 *
            //                   (1e0 + beaconobs->m_values[fi].m_value * 1e-11);

            // Doppler count and delta time (proper)
            const double Ndop =
                beaconobs->m_values[l1i].m_value - pprev_obs->Ls1;
            const auto delta_tau = tproper.delta_sec(pprev_obs->tproper);
            const double NdopDt = Ndop / delta_tau.to_fractional_seconds();
            printf("%.3s Lt(i)=%.6f Lt(i-1)=%.6f", beaconobs->id(),
                   beaconobs->m_values[l1i].m_value, pprev_obs->Ls1);
            printf(" Dt=%.9f", delta_tau.to_fractional_seconds());

            // we will need the true emitter frequency, feT = feN * (1 + Dfe)
            //const double feT =
            //    feN * (1e0 + Filter.rfoff_estimate(receiver_number, tl1));

            // Ionospheric path delay in [m/sec]
            const double Dion = (iers2010::C / feN) *
                                (cDiono - pprev_obs->Diono) /
                                delta_tau.to_fractional_seconds();
            printf(" Dion=%.9f", Dion);

            // Tropospheric delay in [m/sec]
            // get current estimate of Wet Zenith delay
            const double cWzd = Filter.tropo_estimate(receiver_number);
            const double Dtropo = 
                (cDtropo.sum(cWzd) - pprev_obs->Dtropo.sum(cWzd)) /
                delta_tau.to_fractional_seconds();
            printf(" Dtropo=%.9f", Dtropo);

            // Relativistic corrections
            const double Drel = (cDrel - pprev_obs->Drel) /
                                delta_tau.to_fractional_seconds();
            printf(" Drelc=%.9f", Drel);

            // we will need the estimated Δf_e / f_eN
#ifdef RANDOM_RFO
            const double DfefeN = Filter.rfoff_estimate(receiver_number);
#else
            const double DfefeN = Filter.rfoff_estimate(receiver_number, tl1);
#endif
            printf(" Dfe/feN=%.6e feN=%.6f frT=%.6f Rt(i)=%.6f Rt(i-1)=%.6f "
                   "[%.6f %.6f %.6f]\n",
                   DfefeN, feN, frT, rho, pprev_obs->rho(), bxyz_ion(0),
                   bxyz_ion(1), bxyz_ion(2));

            // handle range-rate measurement
            // ---------------------------------------------------------
            const double Dtau = delta_tau.to_fractional_seconds();
            // Warning! Uobs and Utheo have oposite signs, aka Uobs ~= -Utheo
            const double Uobs =
                (iers2010::C / feN) * (feN - frT - NdopDt) + Dion + Drel;
            // theoretical value
            const double Utheo = (1e0 / Dtau) * (rho - pprev_obs->rho()) +
                                 Dtropo -
                                 (iers2010::C / feN) * (NdopDt + frT) * DfefeN;

            const double oc = Uobs + Utheo;
            const double threshold =
                (rstats.count() > 5)
                    ? (3e1 * rstats.stddev() / std::sin(el))
                    : 1e3;
            if (std::abs(oc) < threshold) {

              rstats.update(oc);

              [[maybe_unused]] const Eigen::Matrix<double, 3, 3> R =
                  dso::topocentric_matrix(
                      dso::car2ell<dso::ellipsoid::grs80>(bxyz_ion));

              // partials wrt [x,y,z,Vx,Vy,Vz,Np_1, Np_2, ..., Lwz,Dfe/feN]
              // we actually need the partials of -Utheo (not Utheo), so
              // remember to change signs at the end
              Eigen::VectorXd dHdX = Eigen::VectorXd::Zero(NumParams);
              // dz/dy
              dHdX.block<3, 1>(0, 0) =
                  /*R.transpose() *
                  (r_enu / rho - pprev_obs->s / pprev_obs->rho()) *
                  (1e0 / Dtau);*/
                  (r_enu / rho - pprev_obs->s / pprev_obs->rho()) *
                  (1e0 / Dtau);
              dHdX.block<3, 1>(3, 0) = Eigen::VectorXd::Zero(3);
              // dz/dC_drag
              dHdX(Filter.drag_coef_index()) = 1e0;
#ifdef RANDOM_RFO
              // dz/dDf
              dHdX(Filter.rfoff_index(receiver_number)) =
                  (iers2010::C / feN) * (NdopDt + frT);
#else
              const int _idx = Filter.rfoff_index(receiver_number);
              // dz/dDf -- a term in y = a + b*t
              dHdX(_idx) =
                  - (iers2010::C / feN) * (NdopDt + frT);
              // dz/dDf -- b term
              dHdX(_idx+ 1) = - Filter.deltat(tl1) *
                  (iers2010::C / feN) * (NdopDt + frT);
#endif
              // dz/dLwet
              dHdX(Filter.tropo_index(receiver_number)) =
                  (cDtropo.mfw - pprev_obs->Dtropo.mfw) / Dtau;

              // change signs in the vector
              dHdX *= -1e0;

              const Eigen::VectorXd apriori = Filter.estimates();

              // Filter measurement update
              Filter.observation_update(Uobs, -Utheo, obs_sigma / std::sin(el)*std::sin(el),
                                       dHdX);

              dso::strftime_ymd_hmfs(tl1, dtbuf);
              printf("%s (TAI) %.4s %d %+12.3f %+12.3f %+12.3f "
                     "%+10.6f %+10.6f %+10.6f %+9.6e %9.6f %+9.6f %+9.6f %.3f %.12f "
                     "\n",
                     dtbuf, beacon_it->m_station_id, pprev_obs->arcnr,
                     Filter.estimates()(0), Filter.estimates()(1),
                     Filter.estimates()(2), Filter.estimates()(3),
                     Filter.estimates()(4), Filter.estimates()(5),
#ifdef RANDOM_RFO
                     Filter.rfoff_estimate(receiver_number),
#else
                     Filter.rfoff_estimate(receiver_number, tl1),
#endif
                     Filter.tropo_estimate(receiver_number), Filter.drag_coef(),
                     oc, dso::rad2deg(el), tl1.as_mjd());

              svState.state = Filter.estimates().block<6, 1>(0, 0);

            } else {
              dso::strftime_ymd_hmfs(tl1, dtbuf);
              fprintf(stderr,
                      "WARNING! Skipping observation at %s (TAI) for beacon "
                      "%.4s cause of too big "
                      "O-C value %u/%u : %.3f > %.3f\n",
                      dtbuf, beacon_it->m_station_id, ndop_count_rejected,
                      ndop_count, oc, threshold);
              ++ndop_count_rejected;
            }

            // update previous observation for next Ndop
            pprev_obs->update(tl1, tproper, beaconobs->m_values[l1i].m_value,
                              beaconobs->m_values[l2i].m_value, cDiono, cDtropo,
                              cDrel, r_enu, rho);
            ++ndop_count;

          } // computing Ndop

        } // tropo correction (grid) found 
        else { 
          fprintf(stderr, "Skipping observation, cannot copute tropospheric correction!\n");
        }

        } // elevation > limit

      } // flags ok

      // next beacon observation block for this epoch
      ++beaconobs;
    } // while (beaconobs != it.cblock.end())

    svState.state = Filter.estimates().block<6, 1>(0, 0);

    // transform to keplerian elements, just for printing
    if (1 == 2) {
      Eigen::Matrix<double, 6, 1> inertial_state;
      Eigen::Matrix<double, 3, 3> rc2i, rpom;
      double era, lod;
      assert(!gcrs2itrs(tl1.as_mjd(), Integrator.params->eopLUT, rc2i, era,
                        rpom, lod));
      inertial_state =
          dso::yter2cel(svState.state.block<6, 1>(0, 0), rc2i, era, lod, rpom);
      dso::OrbitalElements kepler = dso::state2elements(GM, inertial_state);
      printf("OE: %s (TAI) %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f\n", dtbuf,
             dso::rad2deg(kepler.semimajor()),
             dso::rad2deg(kepler.eccentricity()),
             dso::rad2deg(kepler.inclination()), dso::rad2deg(kepler.Omega()),
             dso::rad2deg(kepler.omega()), dso::rad2deg(kepler.mean_anomaly()));
    }

    if (tl1.delta_sec(rnx.time_of_first_obs()) >
        dso::cast_to<dso::seconds, dso::nanoseconds>(
            dso::seconds(12 * 60 * 60)))
      break;

  } // for every new data block in the RINEX file

  return 0;
}

int get_rinex_indexes(const dso::DorisObsRinex &rnx, int &l1_idx, int &l2_idx,
                      int &f_idx, int &w1_idx, int &w2_idx) noexcept {
  l1_idx = rnx.get_observation_code_index(
      dso::ObservationCode{dso::ObservationType::phase, 1});

  // index of the 400MHz phase measurement (need for iono-free reduction)
  l2_idx = rnx.get_observation_code_index(
      dso::ObservationCode{dso::ObservationType::phase, 2});

  // index of the F measurement (relative frequency offset)
  f_idx = rnx.get_observation_code_index(
      dso::ObservationCode{dso::ObservationType::frequency_offset});

  // index of the W measurement (power level), frequency 2GHz
  w1_idx = rnx.get_observation_code_index(
      dso::ObservationCode{dso::ObservationType::power_level, 1});

  // index of the W measurement (power level), frequency 400MHz
  w2_idx = rnx.get_observation_code_index(
      dso::ObservationCode{dso::ObservationType::power_level, 2});

  if (f_idx < 0 || l1_idx < 0 || l2_idx < 0 || w1_idx < 0 || w2_idx < 0) {
    fprintf(stderr,
            "[ERROR] Failed to find requested Observation Types in RINEX\'s "
            "observation "
            "types vector! (traceback: %s)\n",
            __func__);
    return 1;
  }

  return 0;
}

Eigen::Matrix<double, 3, 1>
beacon_coordinates(const char *_4charid,
                   const std::vector<dso::BeaconCoordinates> &crdVec) noexcept {
  const auto it = std::find_if(crdVec.cbegin(), crdVec.cend(),
                               [&](const dso::BeaconCoordinates &b) {
                                 return !std::strncmp(b.id, _4charid, 4);
                               });
  assert(it != crdVec.cend());
  double data[3] = {it->x, it->y, it->z};
  return Eigen::Matrix<double, 3, 1>(data);
}

int relativity_corrections(const Eigen::Matrix<double, 6, 1> &sv_state,
                           const Eigen::Matrix<double, 3, 1> &rbeacon,
                           double Re, double GM, double J2, double &Drel_c,
                           double &Drel_r) noexcept {
  // correction for beacon (no J2 term) -- emitter
  const double Bdelta_clock = dso::relativistic_clock_correction(
      rbeacon, Eigen::Matrix<double, 3, 1>::Zero(), GM);

  // correction for satellite (including J2 term) -- receiver
  const double Sdelta_clock = dso::relativistic_clock_correction(
      sv_state.block<3, 1>(0, 0), sv_state.block<3, 1>(3, 0), GM, J2, Re);

  // relativity clock correction, Lemoine 2016
  // TODO why the fuck a 0.5 here?
  // Drel_c = .5e0 * (Sdelta_clock - Bdelta_clock);
  Drel_c = Sdelta_clock - Bdelta_clock;

  // TODO relativity correction for travel time
  Drel_r = 0e0;

  return 0;
}

// TODO this should only be done once, for all stations!
/// @param[in] bxyz_arp Beacon ECEF coordinates at ARP
/// @param[in] beacon The beacon
/// @return ECEF coordinates of the beacon's iono-free phase center
Eigen::Matrix<double, 3, 1>
beacon_arp2ion(const Eigen::Matrix<double, 3, 1> &bxyz_arp,
               const dso::BeaconStation &beacon) noexcept {
  //  transform cartesian to ellipsoidal (antenna RP)
  const Eigen::Matrix<double, 3, 1> lfh =
      dso::car2ell<dso::ellipsoid::grs80>(bxyz_arp);
  const Eigen::Matrix<double, 3, 3> R = dso::topocentric_matrix(lfh);
  return bxyz_arp + R.transpose() * beacon.iono_free_phase_center();
}

int prepare_beacon_coordinates(
    std::vector<dso::BeaconCoordinates> &beaconCrdVec,
    const char *beacon_info_tbl,
    const dso::datetime<dso::nanoseconds> &t) noexcept {

  // load and parse the Beacon Information file (TODO btbl.load_to_memmory can
  // throw)
  dso::BeaconInformationTable btbl(beacon_info_tbl);
  const auto tblvec = btbl.load_to_memmory(t);

  // apply station height to every beacon
  for (auto it = beaconCrdVec.begin(); it != beaconCrdVec.end(); it++) {
    auto bhgt =
        std::find_if(tblvec.cbegin(), tblvec.cend(),
                     [&](const dso::BeaconInformationTableEntry &entry) {
                       return (it->id[0] == entry._4charid[0] &&
                               it->id[1] == entry._4charid[1] &&
                               it->id[2] == entry._4charid[2] &&
                               it->id[3] == entry._4charid[3]);
                     });
    if (bhgt == tblvec.cend()) {
      fprintf(stderr,
              "ERROR. Failed to find entry for beacons %.4s in table file %s\n",
              it->id, beacon_info_tbl);
      return 1;
    }

    // (cartesian) coordinates of beacon in DPOD
    double data[3] = {it->x, it->y, it->z};
    const Eigen::Matrix<double, 3, 1> cartesian(data);
    // (ellipsoidal) coordinates of beacon in DPOD
    const auto lfh = dso::car2ell<dso::ellipsoid::grs80>(cartesian);
    // cartesian-to-topocentric matrix
    const auto R = dso::topocentric_matrix(lfh);
    // ENU eccentricity from tables file
    data[0] = 0e0;
    data[1] = 0e0;
    data[2] = bhgt->_height;
    const Eigen::Matrix<double, 3, 1> denu(data);
    // apply ENU and get cartesian coordinates (now w.r.t beacon RP)
    const Eigen::Matrix<double, 3, 1> arp = cartesian + R.transpose() * denu;
    // swap coordinates in vector
    it->x = arp(0);
    it->y = arp(1);
    it->z = arp(2);
  }

  return 0;
}

int get_tropo(const dso::datetime<dso::nanoseconds> &t,
              const Eigen::Matrix<double, 3, 1> &bxyz, double zd,
              const dso::Gpt3Grid &grid, TropoDetails &Dtrop) noexcept {

  // validate zenith angle
  if (!(zd >= 0e0 && zd <= dso::DPI / 2e0)) {
    fprintf(stderr, "WTF!! Weird zenith angle, is: %.2f\n", dso::rad2deg(zd));
  }
  assert(zd >= 0e0 && zd <= dso::DPI / 2e0);

  // ellipsoidal coordinates of the station; store them in an array
  Eigen::Matrix<double, 3, 1> bell = dso::car2ell<dso::ellipsoid::grs80>(bxyz);
  std::vector<std::array<double, 3>> ellipsoidal(
      1, std::array<double, 3>{bell(0), bell(1), bell(2)});
  // note: ellipsoidal = longitude, latitude, h_ell

  // store GPT3 results here
  std::vector<dso::gpt3_result> g3out;

  // call gpt3_fast to get parameters; store them at g3out (including temporal
  // variations)
  if (dso::gpt3_fast(t, ellipsoidal, 0, grid, g3out)) {
    fprintf(stderr, "[ERROR] Failed to compute gpt3!\n");
    return 20;
  }

  // use VMF3 to compute mapping functions (with height correction)
  double mfh, mfw;
  // if (dso::vmf3(g3out[0].ah, g3out[0].aw, t, bell(1), bell(0), zd, mfh, mfw)) {
  if (dso::vmf3_ht(g3out[0].ah, g3out[0].aw, t, bell(1), bell(0), bell(2), zd,
                   mfh, mfw)) {
    fprintf(stderr, "Failed to compute VMF3\n");
    return 30;
  }

  // use refined saastamnoinen to compute the hydrostatic delay (zenith)
  const double zhd0 = dso::saasthyd(g3out[0].p, bell(1), bell(2));
  // apply VMF3 mapping function to compute hydrostatic delay at given zenith
  // const double Dtropo_hydrostatic = zhd0 * vmf3_res.mfh;

  // use Askne and Nordius to approximate wet dealy in zenith
  const double zwd0 = dso::asknewet(g3out[0].e, g3out[0].Tm, g3out[0].la);
  // apply VMF3 mapping function
  // const double Dtropo_wet = zwd0 * vmf3_res.mfw;

  Dtrop.Lhz = zhd0;
  Dtrop.mfh = mfh;
  Dtrop.Lwz = zwd0;
  Dtrop.mfw = mfw;

  return 0;
}

int get_tropo_vmf(const char *site, const dso::datetime<dso::nanoseconds> &t,
              const Eigen::Matrix<double, 3, 1> &bxyz, double zd,
              dso::SiteVMF3Feed& feed, TropoDetails &Dtrop) noexcept {

  // validate zenith angle
  if (!(zd >= 0e0 && zd <= dso::DPI / 2e0)) {
    fprintf(stderr, "WTF!! Weird zenith angle, is: %.2f\n", dso::rad2deg(zd));
  }
  assert(zd >= 0e0 && zd <= dso::DPI / 2e0);

  // ellipsoidal coordinates of the station; store them in an array
  Eigen::Matrix<double, 3, 1> bell = dso::car2ell<dso::ellipsoid::grs80>(bxyz);
  // note: ellipsoidal = longitude, latitude, h_ell

  // interpolate for site
  dso::vmf3_details::SiteVMF3GRMeteoRecord meteo;
  if (feed.interpolate(site, t, meteo))
    return 1;

  double mfh, mfw;
  if (dso::vmf3(meteo.ah, meteo.aw, t, bell(1), bell(0), dso::deg2rad(zd), mfh,
                mfw)) {
    fprintf(stderr, "ERROR Failed to compute vmf3!\n");
    return 1;
  }

  // use refined saastamnoinen to compute the hydrostatic delay (zenith)
  const double zhd0 = meteo.zhd;
  // apply VMF3 mapping function to compute hydrostatic delay at given zenith
  // const double Dtropo_hydrostatic = zhd0 * vmf3_res.mfh;

  // use Askne and Nordius to approximate wet dealy in zenith
  const double zwd0 = meteo.zwd;
  // apply VMF3 mapping function
  // const double Dtropo_wet = zwd0 * vmf3_res.mfw;

  Dtrop.Lhz = zhd0;
  Dtrop.mfh = mfh;
  Dtrop.Lwz = zwd0;
  Dtrop.mfw = mfw;

  return 0;
}
