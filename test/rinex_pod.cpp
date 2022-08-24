#include "datetime/datetime_write.hpp"
#include "doris_observation_equations.hpp"
#include "doris_rinex.hpp"
#include "doris_utils.hpp"
#include "filters/ekf.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iersc.hpp"
#include "iers2010/tropo.hpp"
#include "integrators.hpp"
#include "sp3/sp3.hpp"
#include "sp3/sv_interpolate.hpp" /* debug mode */
#include "var_utils.hpp"
#include <cstdio>
#include <geodesy/geoconst.hpp>

constexpr const double EleCutOff = 7e0; // elevation cut-off angle, [deg]

// max time difference between two observations to mark a new arc pass
constexpr const long max_sec_for_new_arc = 120L;

// Standard gravitational parameters for Sun and Moon in [km^3 / sec^2]
double GMSun, GMMoon;

// usually using these datetimes ...
using Datetime = dso::datetime<dso::nanoseconds>;

struct TropoDetails {
  double Lhz, mfh;
  double Lwz, mfw;
  double sum() const noexcept { return Lhz * mfh + Lwz * mfw; }
  double sum(double Lwz_) const noexcept { return Lhz * mfh + Lwz_ * mfw; }
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
  int arcnr{0};
  int reinitialize{0}; ///< measurement marked as discontinuity
  SatBeacon(const char *id_, int count_, const Datetime &ttai_, const Datetime &tproper_, double L1,
            double L2, double Diono_, const TropoDetails &Dtropo_, double Drel_,
            const Eigen::Matrix<double, 3, 1> &s_) noexcept
      : count(count_), ttai(ttai_), tproper(tproper_), Ls1(L1), Lu2(L2), Diono(Diono_), Drel(Drel_),
        Dtropo(Dtropo_), s(s_) {
    // std::strncpy(id3c, id_, 3); // fuck the warning!
    std::memcpy(id3c, id_, sizeof(char)*3);
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
    return *this;
  }
  void update(const Datetime &ttai_, const Datetime &tproper_, double L1_, double L2_, double Diono_,
              const TropoDetails &Dtropo_, double Drel_,
              const Eigen::Matrix<double, 3, 1> &s_) noexcept {
    ttai = ttai_;
    tproper = tproper_;
    Ls1 = L1_;
    Lu2 = L2_;
    Diono = Diono_;
    Dtropo = Dtropo_;
    Drel = Drel_;
    s = s_;
  }
  double rho() const noexcept { return s.norm(); }
};

// hold satellite state & time
struct SatelliteState {
  double mjd_tai;
  Eigen::Matrix<double, 6, 1> state;    // state vector at t=ttai in ECEF
  Eigen::Matrix<double, 6, 6> Phi;      // state transition matrix at t=tai

  // ECEF to GCRF
  Eigen::Matrix<double, 6, 1>
  celestial(const dso::EopLookUpTable &elut) const noexcept {
    // Terrestrial to Celestial transformation matrix and derivative for this
    // TAI
    Eigen::Matrix<double, 3, 3> dt2c;
    Eigen::Matrix<double, 3, 3> t2c = dso::itrs2gcrs(mjd_tai, elut, dt2c);

    // transform geocentric state to inertial to propagate orbit
    Eigen::Matrix<double, 6, 1> rcel;
    rcel.block<3, 1>(0, 0) = t2c * state.block<3, 1>(0, 0);
    rcel.block<3, 1>(3, 0) =
        t2c * state.block<3, 1>(3, 0) + dt2c * state.block<3, 1>(0, 0);

    return rcel;
  }

  int integrate(double mjd_target, dso::SGOde &integrator) noexcept {
    // Vector containing state + variational equations size: 6 + 6x6
    // Ref. Frame: inertial
    Eigen::Matrix<double, 6 + 6 * 6, 1> yPhi =
        Eigen::Matrix<double, 6 + 6 * 6, 1>::Zero();
    yPhi.block<6, 1>(0, 0) = celestial(integrator.params->eopLUT);
    {
      int k = 6;
      for (int col = 1; col < 7; col++)
        for (int row = 0; row < 6; row++)
          yPhi(k++) = (col - 1 == row) ? 1e0 : 0e0;
    }

    // t0 for variational equations (TAI)
    integrator.params->mjd_tai = mjd_tai;

    // target t for variational equations; seconds after t0
    double tout = (mjd_target - mjd_tai) * 86400e0;

    // set initial intergation flag
    integrator.flag() = 1;

    // keep solution here (celestial RF at tout)
    Eigen::VectorXd sol(6 + 6 * 6);

    // integrate (in inertial RF), from 0 to tout [sec]
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

    // Terrestrial to Celestial transformation matrix and derivative for this
    // TAI
    Eigen::Matrix<double, 3, 3> dt2c;
    Eigen::Matrix<double, 3, 3> t2c =
        dso::itrs2gcrs(mjd_tai, integrator.params->eopLUT, dt2c);

    // transform inertial to terrestrial
    state.block<3, 1>(0, 0) = t2c.transpose() * sol.block<3, 1>(0, 0);
    state.block<3, 1>(3, 0) = t2c.transpose() * sol.block<3, 1>(3, 0) +
                              dt2c.transpose() * sol.block<3, 1>(0, 0);

    // assign Phi matrix (6x6)
    for (int i = 0; i < 6; i++) {
      Phi.col(i) = yPhi.block<6, 1>(6 * (i + 1), 0);
    }

    return 0;
  }
};

int relativity_corrections(const Eigen::Matrix<double, 6, 1> &sv_state,
                           const Eigen::Matrix<double, 3, 1> &rbeacon,
                           double Re, double GM, double J2, double &Drel_c,
                           double &Drel_r) noexcept;

Eigen::Matrix<double, 3, 1>
beacon_arp2ion(const Eigen::Matrix<double, 3, 1> &bxyz_arp,
               const dso::BeaconStation &beacon) noexcept;

// get the l1, l2 and f indexes off from a RINEX file (instance)
int get_rinex_indexes(const dso::DorisObsRinex &rnx, int &l1, int &l2,
                      int &f, int &w1_idx, int &w2_idx) noexcept;

Eigen::Matrix<double, 3, 1>
beacon_coordinates(const char *_4charid,
                   const std::vector<dso::BeaconCoordinates> &crdVec) noexcept;

int get_tropo(const dso::datetime<dso::nanoseconds> &t,
              const Eigen::Matrix<double, 3, 1> &bxyz, double zd,
              const dso::Gpt3Grid &grid, TropoDetails &Dtrop) noexcept;

Eigen::Matrix<double, 3, 1>
range_rate(const Eigen::Matrix<double, 3, 1> &s_t0,
           const Eigen::Matrix<double, 3, 1> &rsta,
           const Eigen::Matrix<double, 3, 1> &rsat, double dt,
           double &rr_computed, Eigen::Matrix<double, 6, 1> &drrdrv) noexcept;

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
  // Construct the DorisRinex instance rnx
  // -------------------------------------------------------------------------
  if (dso::get_yaml_value_depth2(config, "data", "doris-rinex", buf)) {
    fprintf(stderr, "ERROR. Failed parsing data/rinex file from YAML %s\n",
            argv[1]);
    return 1;
  }
  dso::DorisObsRinex rnx(buf);
#ifdef DEBUG
  rnx.print_metadata();
#endif

  // Initial Orbit
  // -------------------------------------------------------------------------
  // Intial satellite state, get it from the Sp3 using the RINEX's time of
  // first obs
  SatelliteState svState;
  dso::get_yaml_value_depth2(config, "data", "sp3", buf);
  dso::Sp3c sp3(buf);
  dso::Sp3Iterator sp3_iterator(sp3);

  // SP3 Validation Orbit (Should not be needed)
  // -------------------------------------------------------------------------
  // dso::Sp3c sp3(buf);
  // dso::sp3::SatelliteId sv("XXX");
  // sv.set_id(sp3.sattellite_vector()[0].id);
  // dso::SvInterpolator sv_intrp(sv, sp3);

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
    dso::EopFile eopin(buf);
    const int ref_mjd = rnx.ref_datetime().as_mjd();
    const int start = ref_mjd - 4;
    const int end = ref_mjd + 4;
    if (eopin.parse(start, end, eop_lut)) {
      fprintf(stderr, "ERROR. Failed collecting EOP data\n");
      return 1;
    }
  }

  // Gravity
  // -------------------------------------------------------------------------
  // parse degree and order and the requested gravity model into a
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
    error = dso::parse_gravity_model(buf, degree, order, harmonics, true);
  if (error) {
    fprintf(stderr, "ERROR Failed handling gravity field model!\n");
    return 1;
  }

  // Troposphere
  // -------------------------------------------------------------------------
  // read-in the grid file. That is all for now
  dso::get_yaml_value_depth3(config, "troposphere", "gpt3", "grid", buf);
  dso::Gpt3Grid gpt3_grid(buf);

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
  if (extrapolate_sinex_coordinates(buf, rnx.stations(), rnx.ref_datetime(),
                                    beaconCrdVec, true)) {
    fprintf(stderr,
            "ERROR. Failed extracting/extrapolating beacon coordinates\n");
    return 1;
  }

  // Previous Observation relating Beacon/Satellite (to compute Ndop)
  // -------------------------------------------------------------------------
  std::vector<SatBeacon> prevec;
  prevec.reserve(beaconCrdVec.size());

  // Setup Integration Parameters for Orbit Integration
  // We will need the pck (SPICE) kernel for gravitational parameters of Sun
  // and Moon
  // -------------------------------------------------------------------------
  if (dso::get_yaml_value_depth2(config, "naif-kernels", "pck", buf)) {
    fprintf(stderr, "ERROR Failed locating NAIF pck kernel\n");
    return 1;
  }
  dso::IntegrationParameters IntegrationParams(degree, order, eop_lut,
                                               harmonics, buf);

  // Orbit Integrator
  // -------------------------------------------------------------------------
  // Setup an integrator, to extrapolate orbit with:
  // 1. Relative accuracy 1e-12
  // 2. Absolute accuracy 1e-12
  // 3. Num of Equations: 6 for state and 6*6 for variational equations
  dso::SGOde Integrator(dso::VariationalEquations, 6 + 6 * 6, 1e-12, 1e-12,
                        &IntegrationParams);

  // get the (RINEX) indexes for the observables we want
  int l1i, l2i, fi, w1i, w2i;
  if (get_rinex_indexes(rnx, l1i, l2i, fi, w1i, w2i))
    return 1;

  // Extended Kalman Filter instance
  // -------------------------------------------------------------------------
  const int NumParams =
      6                        ///< satellite state
      + rnx.stations().size()  ///< beacon relative frequency offset
      + rnx.stations().size(); ///< wet tropo path dealy (zenith)
  dso::ExtendedKalmanFilter<dso::nanoseconds> Filter(NumParams);

  // Start RINEX data-block iteration
  // -------------------------------------------------------------------------
  // get an iterator to the RINEXs data blocks
  dso::RinexDataBlockIterator it(&rnx);

  // Some variables ...
  const double J2 = harmonics.J2();
  const double GM = harmonics.GM();
  const double Re = harmonics.Re();

  // I only need TLSB at this point; get its RINEX-internal, 3-char id
  const char *tlsb = rnx.beacon_id2internal_id("TLSB");
  assert(tlsb);

  printf("Site Arc Date Time SecOfDay Elev Ndop Dt (c/feN)*(feN-frT-Ndop/dt) (rho2-rho1)/dt Ion Dion Tro Dtro Relc Drelc frT Df/f ClkOffset\n");

  error = 0;
  int dummy_counter = 0;
  // count beacons as we encounter them (above min. elevation)
  int receiver_count = 0;
  // for every new data block in the RINEX file (aka every epoch) ...
  // dso::datetime<dso::nanoseconds> last_ex_target;
  while (!(error = it.next())) {

    // limit to one beacon TLSB; if the block does not cotain this beacon, skip)
    if (it.contains_beacon(tlsb) != it.cblock.end()) {

      // the current reference time for the L1 observation (corrected for
      // receiver clock offset). That is tl1 is approximately TAI.
      // aka proper-time to coordinate-time
      const auto tl1 = it.corrected_l1_epoch();

      // current proper time (aka tau)
      const auto tproper = it.proper_time();

      // a buffer to write datetime strings to
      char dtbuf[64];

      // integrate orbit to here (TAI) TODO
      // svState will contain satellite state for time tl1 in ECEF 
      // first get reference state from sp3, for an epoch as close as possible
      if (sp3_iterator.goto_epoch(tl1)) {
        fprintf(stderr, "ERROR Failed to get reference positio from SP3\n");
        return 1;
      }
      svState.mjd_tai  = sp3_iterator.current_time().as_mjd();
      svState.state(0) = sp3_iterator.data_block().state[0] * 1e3;
      svState.state(1) = sp3_iterator.data_block().state[1] * 1e3;
      svState.state(2) = sp3_iterator.data_block().state[2] * 1e3;
      svState.state(3) = sp3_iterator.data_block().state[4] * 1e-1;
      svState.state(4) = sp3_iterator.data_block().state[5] * 1e-1;
      svState.state(5) = sp3_iterator.data_block().state[6] * 1e-1;
      if (svState.integrate(tl1.as_mjd(), Integrator)) {
      // WARNING just for debugging!
      // Extrapolate orbit to the next second!
      //auto tepoch = sp3_iterator.current_time();
      //tepoch.add_seconds(dso::nanoseconds(1'000'000'000 * 600L));
      //if (tepoch!=last_ex_target) {
      //if (svState.integrate(tepoch.as_mjd(), Integrator)) {
        fprintf(stderr, "ERROR. Failed to integrate orbit!\n");
        return 1;
      }
      //printf("%s %.6f %.6f %.6f %.9f %.9f %.9f %.12f\n",
      //       dso::strftime_ymd_hmfs<dso::nanoseconds>(tepoch, dtbuf),
      //       svState.state(0), svState.state(1), svState.state(2),
      //       svState.state(3), svState.state(4), svState.state(5), tepoch.as_mjd());
      //last_ex_target = tepoch;
      //}

        // if this is the first iteration, initialize the Kalman filter
        if (!dummy_counter) {
          Filter.t = tl1;
          Filter.x.block<6, 1>(0, 0) = svState.state;
          for (int i = 6; i < NumParams; i += 2) {
            Filter.x(i) = 1e-6;    // beacon relative frequency offset
            Filter.x(i + 1) = 0e0; // apriori tropo wet delay at zenith
          }
          // default sigma for releative frequency offset
          Filter.P = Eigen::MatrixXd::Zero(NumParams, NumParams);
          // default sigma for position is 1 [m]
          Filter.P.block<3, 3>(0, 0) =
              5e0 * Eigen::Matrix<double, 3, 3>::Identity();
          // default sigma for velocity is 5 [m/sec]
          Filter.P.block<3, 3>(3, 3) =
              5e0 * Eigen::Matrix<double, 3, 3>::Identity();
        }

        // iterate through the observation set (aka the various beacons with
        // observations for current epoch)
        auto beaconobs = it.cblock.begin();
        while (beaconobs != it.cblock.end()) {

          // an iterator to the RINEX's stations vector (aka a BeaconStation)
          // matching the current beacon
          auto beacon_it =
              rnx.beacon_internal_id2BeaconStation(beaconobs->id());
          assert(beacon_it != rnx.stations().cend());

          // (debugging) only consume TLSB
          if (!std::strncmp(beacon_it->m_station_id, "TLSB", 4)) {

            // check flags
            if (beaconobs->m_values[l1i].m_flag1 ==
                    '1' || // 2 GHz central frequency measurement
                beaconobs->m_values[l1i].m_flag2 ==
                    '1' || // discontinuity of 2 GHZ measurement
                beaconobs->m_values[l2i].m_flag1 ==
                    '1' || // 400 MHz central frequency measurement
                beaconobs->m_values[l2i].m_flag2 ==
                    '1' || // discontinuity of 400 MHZ measurement
                beaconobs->m_values[w1i].m_flag1 ==
                    '1' || // station on restart mode
                beaconobs->m_values[w2i].m_flag1 ==
                    '1' // station on restart mode
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
              // printf("\t\t#Note! Observation has flags:
              // L1[%c%c]/L2[%c%c]/W1[%c%c]/W2[%c%c]; skipped!\n",
              //        beaconobs->m_values[l1i].m_flag1,
              //        beaconobs->m_values[l1i].m_flag2,
              //        beaconobs->m_values[l2i].m_flag1,
              //        beaconobs->m_values[l2i].m_flag2,
              //        beaconobs->m_values[w1i].m_flag1,
              //        beaconobs->m_values[w1i].m_flag2,
              //        beaconobs->m_values[w2i].m_flag1,
              //        beaconobs->m_values[w2i].m_flag2
              //        );

            } else {
              // flags ok, continue with processing ...

              // we are going to need the beacons ECEF coordinates (note that
              // these are antenna RP coordinates)
              const Eigen::Matrix<double, 3, 1> bxyz_sta =
                  beacon_coordinates(beacon_it->m_station_id, beaconCrdVec);

              // Iono-Free phase center w.r.t antenna RP, Cartesian ECEF
              const Eigen::Matrix<double, 3, 1> bxyz_ion =
                  beacon_arp2ion(bxyz_sta, *beacon_it);

              // get azimouth [rad], elevation [rad] and geometric distance [m]
              // (beacon to satellite)
              const Eigen::Matrix<double, 3, 1> r_enu =
                  dso::car2top<dso::ellipsoid::grs80>(
                      bxyz_ion, svState.state.block<3, 1>(0, 0));
              double az, el;
              double rho = dso::top2dae(r_enu, az, el);

              // only process observations to elevation > EleCutOff [deg]
              if (dso::rad2deg(el) < EleCutOff) {
                // dso::strftime_ymd_hmfs(it.cheader.m_epoch, dtbuf);
                // printf("\t\tSkipping observation to beacon %.4s because "
                //        "elevation is "
                //        "%.1f [deg] time: %s (TAI)\n",
                //        beacon_it->m_station_id, dso::rad2deg(el), dtbuf);
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
                      dso::nanoseconds(
                          max_sec_for_new_arc *
                          dso::nanoseconds::sec_factor<
                              dso::nanoseconds::underlying_type>())) {
                    start_new_arc = true;
                  }
                }

                // nominal frequencies for the beacon: feN
                int k; // shift factor
                if (rnx.beacon_shift_factor(beaconobs->id(), k)) {
                  fprintf(stderr,
                          "Failed to find shift factor for beacon %.3s\n",
                          beaconobs->id());
                  return 1;
                }
                double feN, fe2N; // [Hz]
                dso::beacon_nominal_frequency(k, feN, fe2N);

                // ionospheric correction (actually L_2GHz = L_2GHz + Diono)
                const double Diono = dso::carrier_iono_correction(
                    beaconobs->m_values[l1i].m_value,
                    beaconobs->m_values[l2i].m_value);

                // tropospheric correction
                TropoDetails cDtropo;
                if (get_tropo(tl1, bxyz_ion, dso::DPI / 2e0 - el, gpt3_grid,
                              cDtropo)) {
                  return 3;
                }

                // relativistic correction
                double Drel_c, Drel_r;
                relativity_corrections(svState.state, bxyz_ion, Re, GM, J2,
                                       Drel_c, Drel_r);
                const double cDrel = Drel_c + Drel_r;

                // Are we going to compute Doppler count / observations
                // equation? We may not need to, if :
                // 1. this is the first observation of the beacon
                //    aka pprev_obs == prevec.end()
                // 2. We are starting a new arc (aka start_new_arc==1)
                // 3. Previous observation was marked
                // (pprev_obs->reinitialize==1)
                if (pprev_obs == prevec.end() || pprev_obs->reinitialize ||
                    start_new_arc) {
                  if (pprev_obs != prevec.end()) {
                    // not the first observation for the beacon; update last
                    // observation info
                    pprev_obs->update(tl1, tproper,
                                      beaconobs->m_values[l1i].m_value,
                                      beaconobs->m_values[l2i].m_value, Diono,
                                      cDtropo, cDrel, r_enu);
                    // update arc number -- if needed
                    if (start_new_arc)
                      pprev_obs->arcnr += 1;
                    // passed discontinuity, observation ok
                    if (pprev_obs->reinitialize)
                      pprev_obs->reinitialize = 0;
                  } else {
                    // first observation for the beacon
                    prevec.emplace_back(
                        SatBeacon(beaconobs->id(), receiver_count, tl1, tproper,
                                  beaconobs->m_values[l1i].m_value,
                                  beaconobs->m_values[l2i].m_value, Diono,
                                  cDtropo, cDrel, r_enu));
                    ++receiver_count;
                  }

                  // update a-priori value for this beacon zenith wet tropo
                  // delay [m]
                  Filter.x(6 + receiver_count * 2 + 1) = cDtropo.Lwz;

                } else {
                  // compute Ndop and observation equation

                  // we need to find the true proper frequency of the receiver
                  // (aka satellite), f_rT [Hz]
                  const double frT =
                      dso::DORIS_FREQ1_MHZ * 1e6 *
                      (1e0 + beaconobs->m_values[fi].m_value * 1e-11);

                  // Doppler count and delta time (proper)
                  const double Ndop =
                      beaconobs->m_values[l1i].m_value - pprev_obs->Ls1;
                  const auto delta_tau = tproper.delta_sec(pprev_obs->tproper);
                  const double NdopDt =
                      Ndop / delta_tau.to_fractional_seconds();

                  // Ionospheric path delay in [m/sec]
                  const double Dion = (iers2010::C / feN) *
                                      (Diono - pprev_obs->Diono) /
                                      delta_tau.to_fractional_seconds();

                  // Tropospheric delay in [m/sec]
                  // get current estimate of Wet Zenith delay
                  const double cWzd = Filter.x(6 + receiver_count * 2 + 1);
                  [[maybe_unused]] const double Dtropo =
                      (cDtropo.sum(cWzd) - pprev_obs->Dtropo.sum(cWzd)) /
                      delta_tau.to_fractional_seconds();

                  // Relativistic corrections
                  const double Drel = (cDrel - pprev_obs->Drel) /
                                      delta_tau.to_fractional_seconds();

                  dso::strftime_ymd_hmfs(tl1, dtbuf);
                  printf("%.4s %d %s %.6f %.1f %.6f %.6f %.6f %.6f %.6f %.6f "
                         "%.6f %.6f %.6f %.6f %.3f %.6f %.6f\n",
                         beacon_it->m_station_id, pprev_obs->arcnr, dtbuf,
                         tl1.sec().to_fractional_seconds(), dso::rad2deg(el),
                         Ndop, delta_tau.to_fractional_seconds(),
                         (iers2010::C / feN) * (feN - frT - NdopDt),
                         (rho - pprev_obs->rho()) /
                             delta_tau.to_fractional_seconds(),
                         (iers2010::C / feN) * Diono, Dion, cDtropo.sum(),
                         (cDtropo.sum() - pprev_obs->Dtropo.sum()) /
                             delta_tau.to_fractional_seconds(),
                         cDrel, Drel, frT, beaconobs->m_values[fi].m_value,
                         it.cheader.m_clock_offset);

                  // we will need the estimated Δf_e / f_eN
                  const int receiver_number = pprev_obs->count;
                  [[maybe_unused]] const double DfefeN =
                      Filter.x(6 + receiver_number * 2);

                  // Filter time-update
                  // ----------------------------------------------------------------
                  // State transition matrix (augmented)
                  Eigen::MatrixXd PhiP =
                      Eigen::MatrixXd::Identity(NumParams, NumParams);
                  PhiP.block<6, 6>(0, 0) = svState.Phi;
                  /// already / done
                  auto estimates = Filter.x;
                  estimates.block<6, 1>(0, 0) = svState.state;
                  Filter.time_update(tl1, estimates, PhiP);

                  // handle range-rate measurement
                  // ---------------------------------------------------------
                  const double Dtau = delta_tau.to_fractional_seconds();
                  const Eigen::Matrix<double, 3, 3> R = dso::topocentric_matrix(
                      dso::car2ell<dso::ellipsoid::grs80>(bxyz_ion));
                  // observed value
                  const double Uobs =
                      (iers2010::C / feN) * (feN - frT - NdopDt) + Diono + Drel;
                  // theoretical value
                  const double Utheo =
                      (1e0 / Dtau) * (rho - pprev_obs->rho()) +
                      Dtropo - (iers2010::C / feN) * (NdopDt + frT) * DfefeN;
                  // partials wrt [x,y,z,Vx,Vy,Vz,Lwz, Dfe/feN]
                  Eigen::VectorXd dHdX = Eigen::VectorXd::Zero(NumParams);
                  dHdX.block<3, 1>(0, 0) =
                      (1e0 / Dtau) * R.transpose() *
                      ((1e0 / pprev_obs->rho()) * pprev_obs->s -
                       (1e0 / rho) * r_enu);
                  dHdX(6 + receiver_number * 2, 1) =
                      (cDtropo.mfw - pprev_obs->Dtropo.mfw) / Dtau;
                  dHdX(6 + receiver_number * 2 + 1, 1) =
                      -(iers2010::C / feN) * (NdopDt + frT);

                  // Filter measurement update
                  Filter.observation_update(Uobs, Utheo, 5e0 / std::cos(el), dHdX);
          
                  /*
                  Eigen::Matrix<double, 6, 1> drrdrv;
                  double rho_dot; // computed value for range-rate
                  const Eigen::Matrix<double, 3, 1> s = range_rate(
                      pprev_obs->s, bxyz_ion, svState.state.block<3, 1>(0, 0),
                      delta_tau.to_fractional_seconds(), rho_dot, drrdrv);
                  */

                  // update previous observation for next Ndop
                  pprev_obs->update(tl1, tproper,
                                    beaconobs->m_values[l1i].m_value,
                                    beaconobs->m_values[l2i].m_value, Diono,
                                    cDtropo, cDrel, r_enu);

                } // computing Ndop

              } // elevation > limit

            } // if (station == TLSB)

          } // L2GHz flags ok

          // next beacon observation block for this epoch
          ++beaconobs;
        } // for every beacon observation set in epoch

    } // if (it->contains_beacon(tlsb) != it->cblock.end())

    //if (++dummy_counter > 1)
    //   break;
  } // for every new data-block/epoch in the RINEX file

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

/*
int get_satellite_sp3_state(const char *sp3fn, double mjd_tai,
                            Eigen::Matrix<double, 6, 1> &state,
                            double &mjd_state) {
  dso::sp3::SatelliteId sv;
  dso::Sp3DataBlock block, blockp;
  dso::Sp3c sp3(sp3fn);

  if (sp3.num_sats() == 1) {
    sv.set_id(sp3.sattellite_vector()[0].id);
  } else {
    fprintf(stderr, "More than one Satellites in sp3 file. This test program "
                    "is only meant to work with one.\n");
    return 1;
  }

  blockp.t = Datetime::max();

  int error = 0;
  while (!(error = sp3.get_next_data_block(sv, block))) {
    // check the heath status
    if (!block.flag.is_set(dso::Sp3Event::bad_abscent_position)) {
      if (block.t.as_mjd() > mjd_tai) {
        break;
      } else {
        blockp = block;
      }
    }
  }
  if (error)
    return error;

  // accumulate state (m, m/sec), earth-fixed
  state << blockp.state[0] * 1e3, blockp.state[1] * 1e3, blockp.state[2] * 1e3,
      blockp.state[4] * 1e-1, blockp.state[5] * 1e-1, blockp.state[6] * 1e-1;

  // set time of state (TAI)
  mjd_state = blockp.t.as_mjd();

  return 0;
}
*/

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
  // correction for beacon (no J2 term)
  const double Bdelta_clock = dso::relativistic_clock_correction(
      rbeacon, Eigen::Matrix<double, 3, 1>::Zero(), GM);

  // correction for satellite (including J2 term)
  const double Sdelta_clock = dso::relativistic_clock_correction(
      sv_state.block<3, 1>(0, 0), sv_state.block<3, 1>(3, 0), GM, J2, Re);

  //// relativity clock correction, Lemoine 2016
  Drel_c = Sdelta_clock - Bdelta_clock;

  // TODO relativity correction for travel time
  Drel_r = 0e0;

  return 0;
}

Eigen::Matrix<double, 3, 1>
beacon_arp2ion(const Eigen::Matrix<double, 3, 1> &bxyz_arp,
               const dso::BeaconStation &beacon) noexcept {
  // transform cartesian to ellipsoidal (antenna RP)
  auto lfh = dso::car2ell<dso::ellipsoid::grs80>(bxyz_arp);
  // just for debugging, validation
  if (lfh(2) < 0e0) {
    fprintf(stderr, "WTF!! Height is negative\n");
    //printf("Cartesian  : %.3f %.3f %.3f [m]\n", bxyz_arp(0), bxyz_arp(1),
    //       bxyz_arp(2));
    //printf("Ellipsoidal: %+.6f %+.6f %+.3f [deg] and [m]\n",
    //       dso::rad2deg(lfh(0)), dso::rad2deg(lfh(1)), dso::rad2deg(lfh(2)));
    //const auto R = dso::topocentric_matrix(lfh);
    //auto b = bxyz_arp + R.transpose() * beacon.iono_free_phase_center();
    //lfh = dso::car2ell<dso::ellipsoid::grs80>(b);
    //printf("Ellipsoidal: %+.6f %+.6f %+.3f [deg] and [m]\n",
    //       dso::rad2deg(lfh(0)), dso::rad2deg(lfh(1)), dso::rad2deg(lfh(2)));
  }
  const auto R = dso::topocentric_matrix(lfh);
  return bxyz_arp + R.transpose() * beacon.iono_free_phase_center();
}

int get_tropo(const dso::datetime<dso::nanoseconds> &t,
              const Eigen::Matrix<double, 3, 1> &bxyz, double zd,
              const dso::Gpt3Grid &grid, TropoDetails &Dtrop) noexcept {

  // validate zenith angle
  assert(zd >= 0e0 && zd <= dso::DPI / 2e0);

  // ellipsoidal coordinates of the station; store them in an array
  Eigen::Matrix<double, 3, 1> bell = dso::car2ell<dso::ellipsoid::grs80>(bxyz);
  std::vector<std::array<double, 3>> ellipsoidal(
      1, std::array<double, 3>{bell(0), bell(1), bell(2)});

  // store GPT3 results here
  std::vector<dso::gpt3_result> g3out;

  // call gpt3_fast to get parameters; store them at g3out
  if (dso::gpt3_fast(t, ellipsoidal, 1, grid, g3out)) {
    fprintf(stderr, "[ERROR] Failed to compute gpt3!\n");
    return 20;
  }

  // use VMF3 to compute mapping functions
  double mfh, mfw;
  if (dso::vmf3(g3out[0].ah, g3out[0].aw, t, bell(1), bell(0), zd, mfh, mfw)) {
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

/*
Eigen::Matrix<double, 3, 1>
range_rate(const Eigen::Matrix<double, 3, 1> &s_t0,
           const Eigen::Matrix<double, 3, 1> &rsta,
           const Eigen::Matrix<double, 3, 1> &rsat, double dt,
           double &rr_computed, Eigen::Matrix<double, 6, 1> &drrdrv) noexcept {

  // ellipsoidal coordinates at t
  const Eigen::Matrix<double, 3, 1> ell =
      dso::car2ell<dso::ellipsoid::grs80>(rsta);

  // ITRF to topocentric matrix with center at beacon
  const auto E = dso::topocentric_matrix(ell(0), ell(1));

  // satellite-beacon vector in topocentric frame
  const Eigen::Matrix<double, 3, 1> s = E * (rsat - rsta);

  // vellocity in topocentric frame
  const Eigen::Matrix<double, 3, 1> v = (s - s_t0) / dt;

  // norms
  const double sn = s.norm();
  const double vn = v.norm();

  // model value for range-rate
  rr_computed = s.dot(v) / sn;

  // partials w.r.t r (topocentric -> ECEF)
  drrdrv.block<3, 1>(0, 0) =
      (sn * v.transpose() - vn * s.transpose()) * E / sn / sn;

  // partials w.r.t v (topocentric -> ECEF)
  drrdrv.block<3, 1>(3, 0) = s.transpose() * E / sn;

  return s;
}
*/
