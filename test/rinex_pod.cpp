#include "doris_rinex.hpp"
#include "doris_utils.hpp"
#include "var_utils.hpp"
#include <charconv>
#include <cstdio>

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

  // Start RINEX data-block iteration
  // -------------------------------------------------------------------------
  // get an iterator to the RINEXs data blocks
  ids::RinexDataBlockIterator it(&rnx);

  int error = 0;
  // for every new data block in the RINEX file ...
  while (!(error = it.next())) {
    // the current time ...
    auto tnow = it.cheader.m_epoch;
  }

  return 0;
}
