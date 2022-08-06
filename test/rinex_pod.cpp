#include "doris_rinex.hpp"
#include "doris_utils.hpp"
#include "var_utils.hpp"
#include "orbit_integration.hpp"
#include <charconv>
#include <cstdio>
#include <datetime/dtfund.hpp>

// Standard gravitational parameters for Sun and Moon in [km^3 / sec^2]
double GMSun, GMMoon;

// usually using these datetimes ...
using Datetime = dso::datetime<dso::nanoseconds>;

// get the l1, l2 and f indexes off from a RINEX file (instance)
int get_rinex_indexes(const dso::DorisObsRinex &rnx, int &l1, int &l2,
                      int &f) noexcept;


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
  
  // EOP Look Up Table
  // Parse the input EOP data file to create an EopLookUpTable eop_lut
  // -------------------------------------------------------------------------
  dso::EopLookUpTable eop_lut;
  if (dso::get_yaml_value_depth2(config, "eop-info", "eop-file", buf)) {
    fprintf(stderr, "ERROR. Failed parsing eop-info/eop-file file from YAML %s\n",
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

  // station/beacon coordinates
  // Get beacon coordinates from sinex file and extrapolate to RINEX ref. time
  // Result coordinates per beacon are stored in the beaconCrdVec vector.
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
  beaconCrdVec.reserve(70);
  if (extrapolate_sinex_coordinates(buf, rnx.stations(), rnx.ref_datetime(),
                                    beaconCrdVec, true)) {
    fprintf(stderr,
            "ERROR. Failed extracting/extrapolating beacon coordinates\n");
    return 1;
  }

  // Setup Integration Parameters for Orbit Integration
  // ------------------------------------------------------------------------
  dso::IntegrationParameters IntegrationParams(
      degree, order, eop_lut, harmonics);

  // Orbit Integrator
  // -------------------------------------------------------------------------
  // Setup an integrator, to extrapolate orbit with:
  // 1. Relative accuracy 1e-12
  // 2. Absolute accuracy 1e-12
  // 3. Num of Equations: 6 for state and 6*6 for variational equations
  dso::SGOde Integrator(VariationalEquations, 6 + 6*6, 1e-12, 1e-12, &params);


  // get the (RINEX) indexes for the observables we want
  int l1i,l2i,fi;
  if (get_rinex_indexes(rnx, l1i, l2i, fi)) return 1;

  // Start RINEX data-block iteration
  // -------------------------------------------------------------------------
  // get an iterator to the RINEXs data blocks
  dso::RinexDataBlockIterator it(&rnx);

  error = 0;
  // for every new data block in the RINEX file ...
  while (!(error = it.next())) {
    // the current reference time for the L1 observation
    auto tl1 = it.corrected_l1_epoch();
    
    // integrate orbit to current time (TAI) TODO
  }

  return 0;
}

int get_rinex_indexes(const dso::DorisObsRinex &rnx, int &l1_idx, int &l2_idx,
                      int &f_idx) noexcept {
  l1_idx = rnx.get_observation_code_index(
      dso::ObservationCode{dso::ObservationType::phase, 1});

  // index of the 400MHz phase measurement (need for iono-free reduction)
  l2_idx = rnx.get_observation_code_index(
      dso::ObservationCode{dso::ObservationType::phase, 2});

  // index of the F measurement (relative frequency offset)
  f_idx = rnx.get_observation_code_index(
      dso::ObservationCode{dso::ObservationType::frequency_offset});

  if (f_idx < 0 || l1_idx < 0 || l2_idx < 0) {
    fprintf(stderr,
            "[ERROR] Failed to find requested Observation Types in RINEX\'s "
            "observation "
            "types vector! (traceback: %s)\n",
            __func__);
    return 1;
  }

  return 0;
}
