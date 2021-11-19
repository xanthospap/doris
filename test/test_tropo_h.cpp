#include "cppsp3/sv_interpolate.hpp"
#include "doris_rinex.hpp"
#include "doris_utils.hpp"
#include "ggeodesy/car2ell.hpp"
#include "iers2010/tropo.hpp"
#include <algorithm>
#include <cmath>

int beacon_vector2str_array(const std::vector<ids::BeaconStation> &beacons,
                            char **&strarray) noexcept {
  strarray = new char *[beacons.size()];
  for (int i = 0; i < (int)beacons.size(); i++) {
    strarray[i] = new char[4];
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
    std::strncpy(strarray[i], beacons[i].m_station_id, 4);
#pragma GCC diagnostic pop
  }
  return 0;
}

ids::BeaconCoordinates *get_beacon_crd(const char *beacon_4char_id,
                                       ids::BeaconCoordinates *crd_vec,
                                       int num_sta) noexcept {
  for (int i = 0; i < num_sta; i++) {
    if (!std::strncmp(beacon_4char_id, crd_vec[i].id, 4)) {
      return crd_vec + i;
    }
  }
  return nullptr;
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr, "\nUsage: %s [RINEX_FILE] [SINEX] [GPT3_5_GRID]\n", argv[0]);
    return 1;
  }

  // DORIS RINEX instance
  ids::DorisObsRinex rnx(argv[1]);
  // rnx.print_metadata();

  // index of the F measurement (relative frequency offset)
  int p_idx = rnx.get_observation_code_index(
      ids::ObservationCode{ids::ObservationType::ground_pressure});
  if (p_idx < 0) {
    fprintf(stderr,
            "[ERROR] Failed to find ObservationType P in RINEX\'s observation "
            "types vector! (traceback: %s)\n",
            __func__);
    return 1;
  }

  // use the SINEX file to extract beacon coordinates (for reference time of
  // RINEX)
  ids::BeaconCoordinates *crd =
      new ids::BeaconCoordinates[rnx.stations().size()];
  char **sites_array = nullptr;
  beacon_vector2str_array(rnx.stations(), sites_array);
  if (ids::extrapolate_sinex_coordinates(argv[2], sites_array,
                                         rnx.stations().size(),
                                         rnx.ref_datetime(), crd)) {
    fprintf(stderr, "Failed to extrapolate beacon coordinates!\n");
    return 9;
  }
  // dealocate memory
  for (int i = 0; i < (int)rnx.stations().size(); i++)
    delete[] sites_array[i];
  delete[] sites_array;

  // transform beacon coordinates to geodetic, aka [x,y,z]->[lat,lon,hell]
  double lat, lon, hgt;
  for (int i = 0; i < (int)rnx.stations().size(); i++) {
    dso::car2ell<dso::ellipsoid::grs80>(crd[i].x, crd[i].y, crd[i].z, lat, lon,
                                        hgt);
    crd[i].x = lat;
    crd[i].y = lon;
    crd[i].z = hgt;
  }

  // get an iterator to the RINEXs data blocks
  ids::RinexDataBlockIterator it(&rnx);

  int error = 0;
  // get the first data block
  if ((error = it.next())) {
    return error;
  }

    ids::BeaconCoordinates *tls_ptr = nullptr;
  // for every new data block in the RINEX file ...
  while (!(error = it.next())) {
    // the current time ...
    auto tnow = it.cheader.m_epoch;

    auto beaconobs_set = it.cblock.begin(); // Current beacons's observation set
    while (beaconobs_set != it.cblock.end()) {
      // internal 3char name of beacon
      const char *beacon_internal_id = beaconobs_set->m_beacon_id;

      // only process if beacon is Toulousse
      if (!std::strncmp(rnx.beacon_internal_id2id(beacon_internal_id), "TLS",3)) {

      // get the value of pressure for the time/beacon
      double p = beaconobs_set->m_values[p_idx].m_value;
      // compute hydrostatic delay using p observation
      ids::BeaconCoordinates *crd_ptr =
          get_beacon_crd(rnx.beacon_internal_id2id(beacon_internal_id), crd,
                         rnx.stations().size());
      tls_ptr = crd_ptr;
      if (!crd_ptr) {
        fprintf(
            stderr, "[ERROR] Failed to get coordinates for beacon %.3s/%.4s\n",
            beacon_internal_id, rnx.beacon_internal_id2id(beacon_internal_id));
        return 18;
      }
      double hydr_rnx = dso::saasthyd(p, crd_ptr->x, crd_ptr->z);

      // compute hydrostatic delay using gpt3/saastamoinen
      dso::gpt3_result g3out;
      if (dso::gpt3_5_fast(tnow, &crd_ptr->x, &crd_ptr->y, &crd_ptr->z, 1, 0, argv[3], &g3out)) {
          fprintf(stderr, "[ERROR] Failed to compute gpt3!\n");
          return 20;
      }
      double hydr_gpt = dso::saasthyd(g3out.p, crd_ptr->x, crd_ptr->z);

      printf("%.4s %.12f %.5f %.5f %.5f %.5f\n",
             rnx.beacon_internal_id2id(beacon_internal_id), tnow.as_mjd(), p,
             g3out.p, hydr_rnx, hydr_gpt);
      }// beacon is TLS?
      ++beaconobs_set;
    }
  }

printf("Read lines; ended with %d\n", error);
printf("TLS coordinates: lat: %.20f, longitude: %.20f, hell: %.15f\n", tls_ptr->x, tls_ptr->y, tls_ptr->z);

return 0;
}