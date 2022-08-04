#include "doris_rinex.hpp"
#include "doris_utils.hpp"
#include "var_utils.hpp"
#include <charconv>
#include <cstdio>

// approximate number of data points in a bulletin B file (disregarding
// preliminery values)
constexpr const int BBSZ = 40;
// Standard gravitational parameters for Sun and Moon in [km^3 / sec^2]
double GMSun, GMMoon;
// usually using these datetimes ...
using Datetime = dso::datetime<dso::nanoseconds>;
struct EopInfo {
  // actual size of arrays
  int sz;
  // arrays of EOP values extracted from C04
  double mjd[BBSZ], xpa[BBSZ], ypa[BBSZ], ut1a[BBSZ];
};
// to transfer parameters for Variational Equations
struct AuxParams {
  double mjd_tai;
  EopInfo *eopLookUpTables;
  dso::HarmonicCoeffs *hc;
  dso::Mat2D<dso::MatrixStorageType::Trapezoid> *V, *W;
  int degree, order;
};

// get the l1, l2 and f indexes off from a RINEX file (instance)
int get_rinex_indexes(const ids::DorisObsRinex &rnx, int &l1, int &l2,
                      int &f) noexcept;
// get EOPs from a given bulletin-b file for given date. This will fill-up
// the passed-in EopInfo instance (so we can later query it)
int getMeEops(const Datetime &t, char *bulletinb_fn, EopInfo *eopLUT,
              int download = 1) noexcept;


int main(int argc, char *argv[]) {
  // check input
  if (argc != 2) {
    fprintf(stderr, "USAGE: %s [YAML CONFIG]\n", argv[0]);
    return 1;
  }

  // resolve the yaml config file and get the root node
  const YAML::Node config = YAML::LoadFile(argv[1]);
  char buf[256];

  // DORIS RINEX instance
  // -------------------------------------------------------------------------
  if (dso::get_yaml_value_depth2(config, "data", "doris-rinex", buf)) {
    fprintf(stderr, "ERROR. Failed parsing data/rinex filed from YAML %s\n",
            argv[1]);
    return 1;
  }
  ids::DorisObsRinex rnx(buf);
#ifdef DEBUG
  rnx.print_metadata();
#endif

  // Get beacon coordinates from sinex file and extrapolate to RINEX ref. time
  // Reults coordinates per beacon are stored in the beaconCrdVec vector
  // -------------------------------------------------------------------------
  if (dso::get_yaml_value_depth2(config, "reference-frame",
                                 "station-coordinates", buf)) {
    fprintf(stderr,
            "ERROR. Failed parsing reference-frame/station-coordinates filed "
            "from YAML %s\n",
            argv[1]);
    return 1;
  }
  std::vector<ids::BeaconCoordinates> beaconCrdVec;
  beaconCrdVec.reserve(70);
  if (extrapolate_sinex_coordinates(buf, rnx.stations(), rnx.ref_datetime(),
                                    beaconCrdVec, true)) {
    fprintf(stderr,
            "ERROR. Failed extracting/extrapolating beacon coordinates\n");
    return 1;
  }

  // get the (RINEX) indexes for the observables we want
  int l1i,l2i,fi;
  if (get_rinex_indexes(rnx, l1i, l2i, fi)) return 1;

  // Start RINEX data-block iteration
  // -------------------------------------------------------------------------
  // get an iterator to the RINEXs data blocks
  ids::RinexDataBlockIterator it(&rnx);

  int error = 0;
  // for every new data block in the RINEX file ...
  while (!(error = it.next())) {
    // the current reference time for the L1 observation
    auto tl1 = it.corrected_l1_epoch();
    
    // integrate orbit to current time (TAI) TODO
  }

  return 0;
}

int get_rinex_indexes(const ids::DorisObsRinex &rnx, int &l1_idx, int &l2_idx,
                      int &f_idx) noexcept {
  l1_idx = rnx.get_observation_code_index(
      ids::ObservationCode{ids::ObservationType::phase, 1});

  // index of the 400MHz phase measurement (need for iono-free reduction)
  l2_idx = rnx.get_observation_code_index(
      ids::ObservationCode{ids::ObservationType::phase, 2});

  // index of the F measurement (relative frequency offset)
  f_idx = rnx.get_observation_code_index(
      ids::ObservationCode{ids::ObservationType::frequency_offset});

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
