#include "doris_rinex.hpp"
#include "doris_utils.hpp"
#include "cppsp3/sv_interpolate.hpp"
#include <algorithm>

struct LatestBeaconData {
  ids::BeaconObservations obs{};
  dso::datetime<dso::nanoseconds> t;
};

int beacon_vector2str_array(const std::vector<ids::BeaconStation> &beacons,
                            char **&strarray) noexcept {
  strarray =  new char*[beacons.size()];
  for (int i=0; i<(int)beacons.size(); i++) {
    strarray[i] = new char[4];
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
    std::strncpy(strarray[i], beacons[i].m_station_id, 4);
#pragma GCC diagnostic pop
  }
  return 0;
}

int beacon_shift_factor(const char *beacon_id, const ids::DorisObsRinex *rnx,
                        int &k) noexcept {
  //printf("..searching for beacon shift factor .. num stations=%d, capacity=%d\n", (int)rnx->stations().size(), (int)rnx->stations().capacity());
  //std::for_each(rnx->stations().begin(), rnx->stations().end(),[](ids::BeaconStation b){printf("[%.3s]", b.m_internal_code);});
  //printf("\n");

  auto it = std::find_if(
      rnx->stations().begin(), rnx->stations().end(),
      [&](const ids::BeaconStation &b) {
        return !std::strncmp(beacon_id, b.m_internal_code, 3);
      });

  if (it == rnx->stations().end())
    return 1;

  k = it->m_shift_factor;
  return 0;
}

auto match_beacon(const char *beacon_id,
                  std::vector<LatestBeaconData> &prdata) noexcept {
  auto it = std::find_if(prdata.begin(), prdata.end(),
                         [&](const LatestBeaconData &db) {
                           const char *str = db.obs.m_beacon_id;
                           return !std::strncmp(beacon_id, str, 3);
                         });
  return it;
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr, "\nUsage: %s [RINEX_FILE] [SP3 FILE] [SINEX]\n", argv[0]);
    return 1;
  }

  // get an interpolator for the Sp3 file
  dso::Sp3c sp3(argv[2]);
  if (sp3.num_sats() != 1) {
      fprintf(stderr, "More than one satellites found in sp3 file!\n");
      return 5;
  }
  dso::sp3::SatelliteId sv("XXX");
  sv.set_id(sp3.sattellite_vector()[0].id);
  dso::SvInterpolator sv_intrp(sv, sp3);
  double xyz[3] = {0}, dxdydz[3] = {0};

  // DORIS RINEX instance
  ids::DorisObsRinex rnx(argv[1]);
  rnx.print_metadata();

  // latest data records from all beacons
  std::vector<LatestBeaconData> prv_data;
  prv_data.reserve(rnx.stations().size());

  // index of the 2GHz phase measurement
  int l1_idx = rnx.get_observation_code_index(
      ids::ObservationCode{ids::ObservationType::phase, 1});
  if (l1_idx < 0) {
    fprintf(stderr,
            "[ERROR] Failed to find ObservationType L1 in RINEX\'s observation "
            "types vector! (traceback: %s)\n",
            __func__);
    return 1;
  }

  // index of the F measurement (relative frequency offset)
  int f_idx = rnx.get_observation_code_index(
      ids::ObservationCode{ids::ObservationType::frequency_offset});
  if (f_idx < 0) {
    fprintf(stderr,
            "[ERROR] Failed to find ObservationType F in RINEX\'s observation "
            "types vector! (traceback: %s)\n",
            __func__);
    return 1;
  }

  // use the SINEX file to extract beacon coordinates (for reference time of 
  // RINEX)
  //ids::BeaconCoordinates *crd = new ids::BeaconCoordinates[rnx.stations().size()];
  //char **sites_array = nullptr;
  //beacon_vector2str_array(rnx.stations(), sites_array);
  //if (ids::extrapolate_sinex_coordinates(argv[3], sites_array, rnx.stations().size(), rnx.ref_datetime(), crd)) {
  //  fprintf(stderr, "Failed to extrapolate beacon coordinates!\n");
  //  return 9;
  //}

  // get an iterator to the RINEXs data blocks
  ids::RinexDataBlockIterator it(&rnx);

  int error = 0;
  // get the first data block
  if ((error = it.next())) {
    return error;
  }
        
  // for every new data block in the RINEX file ...
  while (!(error = it.next())) {
    // for every beacon observed ...
    double nf2g, nf4m;

    // the current time ...
    auto tnow = it.cheader.m_epoch;

    // satellite position for current epoch
    // TODO need to cast nanoseconds date to microseconds date ....
    dso::datetime<dso::microseconds> now_ms = tnow.cast_to<dso::microseconds>();
    sv_intrp.interpolate_at(now_ms, xyz);

    auto beaconobs_set = it.cblock.begin(); // Cureent beacons's observattion set
    while (beaconobs_set != it.cblock.end()) {
      // internal 3char name of beacon
      const char *beacon_internal_id = beaconobs_set->m_beacon_id;

      // get shift factor and nominal frequency (Hz)
      int k;
      if (beacon_shift_factor(beacon_internal_id, it.rnx, k)) {
        fprintf(stderr,
                "[ERROR] Failed to match beacon shift factor for beacon %.3s\n",
                beacon_internal_id);
        return 101;
      }
      ids::beacon_nominal_frequency(k, nf2g, nf4m);

      // get doppler count (if previous obs exists)
      auto beac_it = match_beacon(beacon_internal_id, prv_data);
      if (beac_it == prv_data.end()) { // no previous obs; just update the array
        prv_data.emplace_back(LatestBeaconData{*beaconobs_set, tnow});
      } else {
        // previous data exists! get Ndop and update array
        double pvalue = beac_it->obs.m_values[l1_idx].m_value;
        double cvalue = beaconobs_set->m_values[l1_idx].m_value;
        double Ndop = cvalue - pvalue;

        // get delta time (receiver time)
        dso::nanoseconds Dt_nsec = tnow.delta_sec(beac_it->t);
        double Dtsec = Dt_nsec.to_fractional_seconds();

        // get the relative frequency offset (Df/f0)
        double Df0 = beaconobs_set->m_values[f_idx].m_value;

        // update the data array
        beac_it->obs = *beaconobs_set;
        beac_it->t = tnow;

        printf("%.3s %.5f %.5f %.5f %.8f %.5f\n", beacon_internal_id, tnow.as_mjd(),
               nf2g, Ndop, Dtsec, Df0);
      }
    ++beaconobs_set;
    }
  }

  printf("Read lines; ended with %d\n", error);
  return 0;
}
