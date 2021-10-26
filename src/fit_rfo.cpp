#include "doris_rinex.hpp"
#include <algorithm>

// Reminder: in the RINEX files, the F value (aka relative frequency offset for
// the receiver's oscillator) is given in units: 1e-11

int get_rinex_rfo(const char *rnx_fn, long &rfo_collected) {
  using namespace ids;

  rfo_collected = 0;
  DorisObsRinex rnx(rnx_fn); // may throw ....
  char line[DorisObsRinex::MAX_RECORD_CHARS];
  RinexDataRecordHeader hdr;
  std::vector<BeaconObservations> obsvec;

  // index of observable F in the list of observation codes for this RINEX
  auto obs_list = rnx.observation_codes();
  auto fit = std::find_if(
      obs_list.begin(), obs_list.end(), [](const ObservationCode &o) {
        return o.type() == ObservationType::frequency_offset;
      });
  if (fit == obs_list.end()) {
    fprintf(stderr,
            "[ERROR] Failed to find observation type F in list of RINEX "
            "observables! (traceback: %s)\n",
            __func__);
    return 1;
  }
  int findex = std::distance(obs_list.begin(), fit);

  int idx = 0;
  while (rnx.stream() &&
         rnx.stream().getline(line, DorisObsRinex::MAX_RECORD_CHARS)) {

    // resolve the data-block header
    if (int status = rnx.resolve_data_epoch(line, hdr); status) {
      fprintf(stderr,
              "[ERROR] Failed to resolve data block header! error=%d "
              "(traceback: %s)\n",
              status, __func__);
      return 1;
    }

    // get the measurements
    if (int status = rnx.read_data_block(hdr, obsvec); status) {
      fprintf(
          stderr,
          "[ERROR] Failed to resolve data block! error=%d (traceback: %s)\n",
          status, __func__);
      return 1;
    }

    // for a single block, all F measurements/values should be the same!
    double rfo = obsvec[0].m_values[findex].m_value;
    for (const auto &bobs : obsvec) {
      if (bobs.m_values[findex].m_value != rfo) {
        fprintf(stderr,
                "[ERROR] F value differs within the same data block! "
                "(traceback: %s)\n",
                __func__);
        return 1;
      }
    }

    // push back the reolved time and offest
    // rfo[offset + idx] = rfo;
    // t[offset + idx] = hdr.m_epoch;
    printf("%.8f %.8f\n", hdr.m_epoch.as_mjd(), rfo);
    ++idx;
  }

  rfo_collected = idx;

  if (rnx.stream().eof()) {
    rnx.stream().clear();
    return 0;
  }

  return 1;
}

int ids::fit_relative_frequency_offset(char **rinex_fns, int num_rinex) noexcept {

  long rfo_collected = 0;
  int error = 0;

  // for every rinex in the list ...
  for (int i = 0; i < num_rinex; i++) {
    try {
      error = get_rinex_rfo(rinex_fns[i], rfo_collected);
      if (error) return error;
    } catch (std::exception &e) {
      fprintf(stderr,
              "[ERROR] Failed to collectd/update F values for RINEX %s "
              "(traceback: %s)\n",
              rinex_fns[i], __func__);
      return 1;
    }
  }

  return 0;
}