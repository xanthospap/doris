#include "doris_rinex.hpp"
#include <algorithm>

int get_rinex_rfo(const char *rnx_fn) {
    using namespace ids;

    DorisObsRinex rnx(rnx_fn); // may throw ....
    char line[DorisObsRinex::MAX_RECORD_CHARS];
    RinexDataRecordHeader hdr;
    std::vector<BeaconObservations> obsvec;
    dso::datetime<dso::nanoseconds> last_epoch;

    // index of observable F in the list of observation codes for this RINEX
    auto obs_list = rnx.observation_codes();
    auto findex = std::find_if(obs_list.begin(), obs_list.end(), [](const ObservationCode& o){ o.type() == ObservationType::frequency_offset)});
    if (findex == obs_list.end()) {
        fprintf(stderr, "[ERROR] Failed to find observation type F in list of RINEX observables! (traceback: %s)\n", __func__);
        return 1;
    }

  while (rnx.stream() && rnx.stream().getline(line, DorisObsRinex::MAX_RECORD_CHARS)) {
    
    if (rnx.stream().eof())
      break;

    // resolve the data-block header    
    if (int status = rnx.resolve_data_epoch(line, hdr); status) {
        fprintf(stderr, "[ERROR] Failed to resolve data block header! error=%d (traceback: %s)\n", status, __func__);
      return 1;
    }

    // get the measurements
    if (int status = rnx.read_data_block(hdr, obsvec); status) {
        fprintf(stderr, "[ERROR] Failed to resolve data block! error=%d (traceback: %s)\n", status, __func__);
      return 1;
    }

    // for a single block, all F measurements/values should be the same!
    double rfo = obsvec[0].m_values[findex];
    for (const auto& bobs : obsvec) {
        if (bobs.m_values[findex] != rfo) {
            fprintf(stderr, "[ERROR] F value differs within the same data block! (traceback: %s)\n", __func__);
            return 1;
        }
    }
}

int fit_relative_frequency_offset(char **rinex_fn, int num_rinex) noexcept {

    // for every rinex in the list ...
    for (int i=0; i<num_rinex; i++) {
        // create a rinex instance; open file
        ids::DorisObsRinex rnx(rinex_fn[i]);

        // parse the field 'F' for every epoch

    }
}