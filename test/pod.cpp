#include "astrodynamics.hpp"
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
#include "iers2010/cel2ter.hpp"
#include "iers2010/hardisp.hpp"
#include "iers2010/iers2010.hpp"
#include "iers2010/iersc.hpp"
#include "iers2010/tropo.hpp"
#include "integrators.hpp"
#include "satellites/jason3.hpp"
#include "satellites/jason3_quaternions.hpp"
#include "sp3/sp3.hpp"
#include "sp3/sv_interpolate.hpp" /* debug mode */
#include "tides.hpp"
#include "var_utils.hpp"
#include <cassert>
#include <cstdio>
#include <cstring>
#include <datetime/dtcalendar.hpp>
#include <datetime/dtfund.hpp>
#include <iers2010/eop.hpp>

constexpr const int m = 0;
constexpr const int n = 6 + m;

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

// hold satellite state & time
struct OrbitIntegrator {
  // Current datetime in TAI
  dso::TwoPartDate mjd_tai;
  // state vector at t=tai in ECEF, at DORIS receiver RP (Iono-Free)
  Eigen::Matrix<double, 6, 1> state;
  // state transition matrix at t=tai
  Eigen::Matrix<double, 6, 6> Phi;

  Eigen::Matrix<double, 6, 6> extractStateTransitionMatrix(
      const Eigen::Matrix<double, (6 + 6 * n), 1> &yP) noexcept {
    Eigen::Matrix<double, 6, 6> F;
    int of = 6;
    F.row(0) = yP.block<6,1>(of,0);
    F.row(1) = yP.block<6,1>(of+1*n,0);
    F.row(2) = yP.block<6,1>(of+2*n,0);
    F.row(3) = yP.block<6,1>(of+3*n,0);
    F.row(4) = yP.block<6,1>(of+4*n,0);
    F.row(5) = yP.block<6,1>(of+5*n,0);
    return F;
  }

  void setInitialconditions(Eigen::Matrix<double,(6+6*n),1> &yP0) noexcept {
    int of = 6;
    yP0.block<6*n,1>(6,0) = Eigen::Matrix<double,(6*n),1>::Zero();
    for (int i=0; i<6; i++) {
      int start = of + i*n;
      yP0(start+i) = 1e0;
    }
  }

  int integrate(const dso::TwoPartDate &mjd_target,
                dso::SGOde &integrator) noexcept {
    // count calls
    static int call_nr = 0;

    // number of equations: 
    // 6 for the state + 6*n for the variational equations, where n = 6 + m
    constexpr const int NumEqn = 6 + 6*n;

    // an instance to transform between ITRF and GCRF coordinates
    dso::Itrs2Gcrs Rot(mjd_tai.tai2tt().normalized(),
                       &(integrator.params->eop_lookup_table()));

    // Vector containing state + variational equations size: 6 + 6xn
    Eigen::Matrix<double, NumEqn, 1> yPhi0 =
        Eigen::Matrix<double, NumEqn, 1>::Zero();

    // transform state from ITRF to GCRF
    yPhi0.block<6, 1>(0, 0) = Rot.itrf2gcrf(state);

    // Initial condition for state transition matrix Φ(t0,t0) = I
    setInitialconditions(yPhi0);

    // t0 for integration (TAI), aka current reference time
    integrator.params->reference_epoch() = mjd_tai.normalized();

    // target t for integration: seconds after t0
    double tout =
        mjd_target.diff<dso::DateTimeDifferenceType::FractionalSeconds>(
            integrator.params->reference_epoch());

    // integration solution 
    Eigen::VectorXd sol(NumEqn);
  
    // integrate (in inertial RF), from 0+mjd_tai to tout+mjd_tai [sec]
    double tsec = 0e0;
    integrator.flag() = dso::SGOde::IFLAG::RESTART;
    if (integrator.de(tsec, tout, yPhi0, sol) != dso::SGOde::IFLAG::SUCCESS) {
      fprintf(stderr, "[ERROR] Integrator error!\n");
      return 1;
    }

    // after integration, we reached the epoch: t0+tsec[sec] (TAI)
    const dso::TwoPartDate tout_mjd(
        dso::TwoPartDate(mjd_tai + dso::TwoPartDate(0e0, tsec / 86400e0))
            .normalized());

    // let's see were we are at
    if (std::abs(tout_mjd.diff<dso::DateTimeDifferenceType::FractionalSeconds>(
            mjd_target)) > 1e-9) {
      fprintf(stderr,
              "ERROR wanted integration to %.9f and got up to %.9f, that is "
              "%.9f sec apart!\n",
              mjd_target.mjd(), tout_mjd.mjd(),
              (mjd_target.mjd() - mjd_tai.mjd()) * 86400e0);
      return 1;
    }

    // everything seems ok, update state and time
    mjd_tai = tout_mjd;

    // Transform state back to ITRF (from GCRF)
    Rot.prepare(mjd_tai.tai2tt());
    state = Rot.gcrf2itrf(Eigen::Matrix<double,6,1>(sol.block<6,1>(0,0)));

    // assign Phi matrix 6x6
    Phi = extractStateTransitionMatrix(sol);

    ++call_nr;
    return 0;
  }
};

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
  if (dso::get_yaml_value_depth2<double>(config, "filtering",
                                         "elevation-cut-off", EleCutOff)) {
    fprintf(stderr, "ERROR failed to find spk kernel\n");
    return 1;
    printf("[EDBUG] Using elevation cut-off angle = %.3f\n", EleCutOff);
  }

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

  // Construct the Doris RINEX instance rnx
  dso::DorisObsRinex rnx(buf);

  // Initial Orbit
  // -------------------------------------------------------------------------
  // Initial satellite state, get it from the Sp3 using the RINEX's time of
  // first obs
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
    const int start = ref_mjd - 8;
    const int end = ref_mjd + 10;
    // parse C04 EOPs 
    if (parse_iers_C0420(buf, start, end, eop_lut)) {
      fprintf(stderr, "ERROR. Failed collecting EOP data\n");
      return 1;
    }
    // transform the tabulated values (time) to the TT timescale
    eop_lut.utc2tt();
    // regularize ΔUT1 and LOD values (to ΔUT1-R and LOD-R)
    eop_lut.regularize(false);
  }

  // Gravity
  // -------------------------------------------------------------------------
  // parse degree, order and the requested gravity model into a
  // HarmonicCoeffs instance. Note that to compute potential we will need
  // Lagrange polynomials (later on)
  int degree, order;
  {
    error = 0;
    error =
        dso::get_yaml_value_depth2<int>(config, "gravity", "degree", degree);
    error += dso::get_yaml_value_depth2<int>(config, "gravity", "order", order);
  }
  dso::HarmonicCoeffs harmonics(degree);
  {
    if (!error)
      error = dso::get_yaml_value_depth2(config, "gravity", "model", buf);
    // parse the un-normalized harmonic ceofficients from model (gfc format)
    if (!error)
      error = dso::parse_gravity_model(buf, degree, order, rnx.ref_datetime(),
                                       harmonics, false);
    if (error) {
      fprintf(stderr, "ERROR Failed handling gravity field model!\n");
      return 1;
    }
  }

  // Station/Beacon coordinates
  // -------------------------------------------------------------------------
  // Get beacon coordinates from sinex file and extrapolate to RINEX ref. time
  // Result coordinates per beacon are stored in the beaconCrdVec vector.
  // Note that these posotion vectors are w.r.t the beacon/antenna reference
  // point. When in actual processing, this has to be changed, if we are
  // considering iono-free analysis
  std::vector<dso::BeaconCoordinates> beaconCrdVec;
  beaconCrdVec.reserve(rnx.stations().size());
  {
    if (dso::get_yaml_value_depth2(config, "reference-frame",
                                   "station-coordinates", buf)) {
      fprintf(stderr,
              "ERROR. Failed parsing reference-frame/station-coordinates file "
              "from YAML %s\n",
              argv[1]);
      return 1;
    }
    // coordinates of beacons (on antenna RP), ECEF/pdop
    if (extrapolate_sinex_coordinates(buf, rnx.stations(), rnx.ref_datetime(),
                                      beaconCrdVec, false, true)) {
      fprintf(stderr,
              "ERROR. Failed extracting/extrapolating beacon coordinates\n");
      return 1;
    }
    printf("[DEBUG] Number of beacons in RINEX file %d, extracted coordinates "
           "for %d\n",
           (int)rnx.stations().size(), (int)beaconCrdVec.size());
  }

  // Ocean Tides Deformation
  dso::get_yaml_value_depth2(config, "ocean-tides", "blq", buf);
  std::vector<dso::BlqSiteInfo> blqInfoVec;
  dso::Hardisp ocdeform;
  {
    // a vector containing all 4-char site names, to extract BLQ info for
    char *namepool = new char[rnx.stations().size()*5];
    std::memset(namepool, '\0', rnx.stations().size()*5);
    int i = 0;
    for (const auto &s : rnx.stations()) {
      std::memcpy(namepool+i*5, s.m_station_id, sizeof(char)*4);
      ++i;
    }
    std::vector<const char *> sites;
    sites.reserve(rnx.stations().size());
    for (int j=0; j<(int)rnx.stations().size(); j++)
      sites.push_back(namepool+j*5);
    if (dso::parse_blq(buf, blqInfoVec, &sites)) {
      fprintf(stderr, "Failed reading BLQ file %s\n", buf);
      return 1;
    }
    delete[] namepool;
  }

  // Ocean Tides Geopotential
  int oc_degree, oc_order;
  std::vector<dso::DoodsonOceanTideConstituent> vdds;
  dso::get_yaml_value_depth2(config, "ocean-tides", "harmonics", buf);
  error = dso::get_yaml_value_depth2<int>(config, "ocean-tides", "degree", oc_degree);
  error += dso::get_yaml_value_depth2<int>(config, "ocean-tides", "order", oc_order);
  if (error || (oc_order > oc_degree)) {
    fprintf(stderr, "Invalid ocean tide degree/order, %d/%d!\n", oc_degree, oc_order);
    return 1;
  }
  if (dso::memmap_octide_coefficients(buf, vdds, oc_degree, oc_order, 3, 1e-11)) {
    fprintf(stderr, "Failed reading ocean tide loading file %s!\n", buf);
    return 1;
  }
  // An ocean tide instance
  dso::OceanTide octide(vdds, harmonics.GM(), harmonics.Re(), oc_degree, oc_order);

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

  // Orbit Integrator
  // -------------------------------------------------------------------------
  // Setup an integrator, to extrapolate orbit with:
  // 1. Relative accuracy 1e-12
  // 2. Absolute accuracy 1e-12
  // 3. Num of Equations: 6 for state and 6*6 for variational equations
  dso::SGOde Integrator(dso::VariationalEquations2, (6+6*n), 1e-12,
                        1e-15, &IntegrationParams);

  // get the (RINEX) indexes for the observables we want
  int l1i, l2i, fi, w1i, w2i;
  {
  if (get_rinex_indexes(rnx, l1i, l2i, fi, w1i, w2i))
    return 1;
  }

  // Important !! OC-TIDES
  IntegrationParams.octide = &octide;
  // Setup Solid Earth Tide
  dso::SolidEarthTide setide(harmonics.GM(), harmonics.Re(),
                             IntegrationParams.GMMon * 1e9,
                             IntegrationParams.GMSun * 1e9);
  IntegrationParams.setide = &setide;

  OrbitIntegrator svState;

  // Start RINEX data-block iteration
  // -------------------------------------------------------------------------
  // get an iterator to the RINEXs data blocks
  dso::RinexDataBlockIterator it(&rnx);

  // Some variables ...
  [[maybe_unused]]const double J2 = harmonics.J2();
  [[maybe_unused]]const double GM = harmonics.GM();
  [[maybe_unused]]const double Re = harmonics.Re();
  // form a rotation instance (ITRF-to-GCRF) for current epoch
  dso::Itrs2Gcrs Rot(rnx.time_of_first_obs(), &eop_lut);
  // a buffer to write datetime strings to ...
  char dtbuf[64];
  // counters ...  
  [[maybe_unused]]unsigned flaged_obs = 0; // number of observations with 'bad' flags
  [[maybe_unused]]unsigned num_obs = 0;    // observation count, regardless if usable or not
  [[maybe_unused]]unsigned num_blocks = 0; // RINEX block count
  // error flag
  error = 0;
  
  // for every new data block in the RINEX file (aka every epoch) ...
  while (!(error = it.next())) {
    // current proper time (aka τ)
    [[maybe_unused]]const auto tobs_proper    = dso::TwoPartDate(it.proper_time());
    [[maybe_unused]]const auto tobs_proper_dt = it.proper_time();

    // get current observation epoch (tobs) from proper time to TAI
    [[maybe_unused]]const auto tobs_tai_dt  = it.corrected_l1_epoch();
    [[maybe_unused]]const auto tobs_tai = dso::TwoPartDate(tobs_tai_dt);
    
    // get current observation epoch in UTC
    [[maybe_unused]]const auto tobs_utc = tobs_tai.tai2utc();

    dso::strftime_ymd_hmfs(tobs_tai_dt, dtbuf);
    printf("# New observation (RINEX) at %s\n", dtbuf);

    // integrate orbit to here (TAI)
    if (!num_blocks) {
      // for first iteration, get reference state from sp3, for an epoch as 
      // close as possible
      if (sp3_iterator.goto_epoch(tobs_tai_dt)) {
        fprintf(stderr, "ERROR Failed to get reference position from SP3\n");
        return 1;
      }
      svState.mjd_tai =  dso::TwoPartDate(sp3_iterator.current_time());
      svState.state(0) = sp3_iterator.data_block().state[0] * 1e3;
      svState.state(1) = sp3_iterator.data_block().state[1] * 1e3;
      svState.state(2) = sp3_iterator.data_block().state[2] * 1e3;
      svState.state(3) = sp3_iterator.data_block().state[4] * 1e-1;
      svState.state(4) = sp3_iterator.data_block().state[5] * 1e-1;
      svState.state(5) = sp3_iterator.data_block().state[6] * 1e-1;
      if (svState.integrate(tobs_tai, Integrator)) {
        fprintf(stderr, "ERROR. Failed to integrate orbit!\n");
        return 1;
      }
    } else {
      if (svState.integrate(tobs_tai, Integrator)) {
        fprintf(stderr, "ERROR. Failed to integrate orbit!\n");
        return 1;
      }
    }
    dso::strftime_ymd_hmfs(tobs_tai_dt, buf);
    printf("%s %+.6f %+.6f %+.6f %+.9f %+.9f %+.9f\n", buf,
           svState.state(0), svState.state(1), svState.state(2),
           svState.state(3), svState.state(4), svState.state(5));
    ++num_blocks; /* augment data block counter */
  } /* data blocks ended, no more data in RINEX */

  printf("#[DEBUG] Number of data block read: %u\n", num_blocks);
  return 0;
}
