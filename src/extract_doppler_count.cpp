#include "doris_rinex.hpp"
#include "doris_system_info.hpp"
#include <cstring>
#include <algorithm>
#include <ggdatetime/dtfund.hpp>
#ifdef DEBUG
#include "ggdatetime/datetime_write.hpp"
#endif

struct BeaconObservationData {
  ids::BeaconObservations obs;
  dso::datetime<dso::nanoseconds> t;
};

int compute_ndop(const std::vector<ids::BeaconObservations> &obsvec,
    const dso::datetime<dso::nanoseconds> &t,
    std::vector<BeaconObservationData> &data) noexcept {
  const char *newid = nullptr;
  for (const auto &newobs : obsvec) {
    newid = newobs.m_beacon_id;

    auto it = std::find_if(data.begin(), data.end(),
        [&](const BeaconObservationData &bod) noexcept {
        return !std::strcmp(bod.obs.m_beacon_id, newid);
        });
    if (it == data.end()) return 1;

    double dl1 = newobs.m_values[0].m_value - it->obs.m_values[0].m_value;
    double dl2 = newobs.m_values[1].m_value - it->obs.m_values[1].m_value;
    // double dsec = (t.as_undelying_type() - it->t.as_underlying_type()) * dso::nanoseconds::sec_factor<double>();
    auto tdsec = t.delta_sec(it->t);
    double dsec = static_cast<double>(tdsec.as_underlying_type()) / dso::nanoseconds::sec_factor<double>();

    printf("%s %.5f %.5f %.7f\n", newid, dl1, dl2, dsec);
  }

  return 0;
}


int update(const std::vector<ids::BeaconObservations> &obsvec,
           const dso::datetime<dso::nanoseconds> &t,
           std::vector<BeaconObservationData> &data) noexcept {
  const char *newid = nullptr;
  for (const auto &newobs : obsvec) {
    // do we already have a record for this beacon?
    newid = newobs.m_beacon_id;
    auto it = std::find_if(data.begin(), data.end(),
                           [&](const BeaconObservationData &bod) noexcept {
                             return !std::strcmp(bod.obs.m_beacon_id, newid);
                           });
    // if not, append
    if (it == data.end()) {
      data.emplace_back(BeaconObservationData{newobs, t});
    // else set current observation records for this beacon
    } else {
      it->obs = newobs;
      it->t = t;
    }
  }

  return 0;
}

int ids::DorisObsRinex::get_doppler_counts() noexcept {

  long parsed_blocks = 0;

  /* find indexes of L1 and L2 in this RINEX's observation codes vector
  auto l1_ptr = std::find_if(m_obs_codes.cbegin(), m_obs_codes.cend(), ObservationType{ObservationType::phase, 1});
  if (l1_ptr == m_obs_codes.cend()) {
    fprintf(stderr, "[ERROR] Failed to find ObservationType L1 in RINEX\'s observation types vector! (traceback: %s)\n", __func__);
    return 1;
  }
  auto l2_ptr = std::find_if(m_obs_codes.cbegin(), m_obs_codes.cend(), ObservationType{ObservationType::phase, 2});
  if (l2_ptr == m_obs_codes.cend()) {
    fprintf(stderr, "[ERROR] Failed to find ObservationType L2 in RINEX\'s observation types vector! (traceback: %s)\n", __func__);
    return 1;
  }

  // index of L1 and L2 in the m_obs_codes vector
  int l1_idx = std::distance(m_obs_codes.cbegin(), l1_ptr);
  int l2_idx = std::distance(m_obs_codes.cbegin(), l2_ptr);
  */

  // vector of beacon observations; last observation extracted for each beacon
  std::vector<BeaconObservationData> latest_data;
  latest_data.reserve(m_stations.size());

  m_stream.seekg(m_end_of_head);
  char line[MAX_RECORD_CHARS];
  RinexDataRecordHeader hdr;
  std::vector<BeaconObservations> obsvec;

  // read first line ....
  if (!m_stream.getline(line, MAX_RECORD_CHARS)) {
    fprintf(stderr, "Failed to read lines from RINEX file (traceback: %s)\n", __func__);
    return 1;
  }

  // this should be a data record line
  if (int error = resolve_data_epoch(line, hdr); error) {
    fprintf(stderr, "Failed parsing data header line (#1)! (traceback: %s)\n", __func__);
    fprintf(stderr, "Line is: \"%s\"\n", line);
    return 2;
  }

  // now this should be the data
  if (int error=read_data_block(hdr, obsvec); error) {
    fprintf(stderr, "Failed parsing data block! (#2) (traceback: %s)\n", __func__);
    return 3;
  }

  // update latest data for beacons
  update(obsvec, hdr.m_epoch, latest_data);

  // l1_ref = ref_obsvec.m_values[l1_idx].m_value;
  // l2_ref = ref_obsvec.m_values[l2_idx].m_value;

  while (m_stream.getline(line, MAX_RECORD_CHARS)) {
    
    // resolve data block header ...
    if (int error = resolve_data_epoch(line, hdr); error) {
      fprintf(stderr, "Failed parsing data header line (#1)! (traceback: %s)\n", __func__);
      fprintf(stderr, "Line is: \"%s\"\n", line);
      return 2;
    }
    // read the corresponding data block ...
    if (int error=read_data_block(hdr, obsvec); error) {
      fprintf(stderr, "Failed parsing data block! (#2) (traceback: %s)\n", __func__);
      return 3;
    }
  
    // l1_nxt = nxt_obsvec.m_values[l1_idx].m_value;
    // l2_nxt = nxt_obsvec.m_values[l2_idx].m_value;

    // // compute the Doppler count in L1 and L2
    // double ndpo_l1 = next_obsvec.m_values[l1_idx].m_value - l1_ref;
    // double ndpo_l2 = next_obsvec.m_values[l2_idx].m_value - l2_ref;

    // compute Doppler ...
    if (compute_ndop(obsvec, hdr.m_epoch, latest_data)) {
      fprintf(stderr,"[ERROR] Failed computing ndop (traceback: %s)\n", __func__);
      return 5;
    }

    // update latest data
    update(obsvec, hdr.m_epoch, latest_data);
    ++parsed_blocks;
  }
  
  // should have reached EOF ...
  if (m_stream.eof()) {
    printf("Number of data blocks parsed: %ld; EOF reached\n", parsed_blocks);
    return 0;
  }

  fprintf(stderr, "[ERROR] Failed getting new line while not EOF encountered (traceback: %s)\n", __func__);
  return 1;
}
