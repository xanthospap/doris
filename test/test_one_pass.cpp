#include "antenna_pcv.hpp"
#include "doris_rinex.hpp"
#include "doris_utils.hpp"
#include "geodesy/car2ell.hpp"
#include "geodesy/car2top.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "iers2010/tropo.hpp"
#include "sp3/sv_interpolate.hpp"
#include <algorithm>
#include <cmath>

int apply_iono_correction = 1;
int new_arc_on_discontinuity = 1;

struct LatestBeaconData {
  ids::BeaconObservations obs {};
  dso::datetime<dso::nanoseconds> t;
  double sv_beacon_rho; // distance between satellite and beacon at time t
  double l2gh_iono_cor; // iono correction
  int start_new_arc { 0 };
};

void update_beacon_obs(
    const char* beacon_id, std::vector<std::pair<long, long>>& beacon_obs,
    const std::vector<ids::BeaconStation>& beacons) noexcept
{
  auto it = std::find_if(beacons.begin(), beacons.end(),
      [&](const ids::BeaconStation& sta) {
        return !std::strncmp(sta.m_station_id, beacon_id, 4);
      });
  assert(it != beacons.end());
  int idx = std::distance(beacons.begin(), it);
  beacon_obs[idx].first++;
  return;
}
void update_missed_beacon_obs(
    const char* beacon_id, std::vector<std::pair<long, long>>& beacon_obs,
    const std::vector<ids::BeaconStation>& beacons) noexcept
{
  auto it = std::find_if(beacons.begin(), beacons.end(),
      [&](const ids::BeaconStation& sta) {
        return !std::strncmp(sta.m_station_id, beacon_id, 4);
      });
  assert(it != beacons.end());
  int idx = std::distance(beacons.begin(), it);
  beacon_obs[idx].second++;
  return;
}

ids::BeaconCoordinates* get_beacon_crd(const char* beacon_4char_id,
    ids::BeaconCoordinates* crd_vec,
    int num_sta) noexcept
{
  for (int i = 0; i < num_sta; i++) {
    if (!std::strncmp(beacon_4char_id, crd_vec[i].id, 4)) {
      return crd_vec + i;
    }
  }
  return nullptr;
}

int beacon_sv_distance(const char* beacon_4charid,
    const ids::BeaconCoordinates* bpos, int num_beacons,
    const double* svpos, double& distance) noexcept
{
  // auto it =
  //    std::find_if(bpos.begin(), bpos.end(), [&](const ids::BeaconCoordinates
  //    &pos) {
  //      return !(std::strncmp(beacon_4charid, pos.id, 4));
  //    });
  // if (it == bpos.end()) {
  //  fprintf(stderr, "[ERROR] Failed to find beacon coordinates for %.4s\n",
  //          beacon_4charid);
  //  return 1;
  //}
  int index = -1;
  for (int i = 0; i < num_beacons; i++) {
    if (!(std::strncmp(beacon_4charid, bpos[i].id, 4))) {
      index = i;
      break;
    }
  }

  if (index == -1) {
    fprintf(stderr, "[ERROR] Failed to find beacon coordinates for %.4s\n",
        beacon_4charid);
    return 1;
  }

  double dx = svpos[0] - bpos[index].x * 1e3;
  double dy = svpos[1] - bpos[index].y * 1e3;
  double dz = svpos[2] - bpos[index].z * 1e3;
  distance = std::sqrt(dx * dx + dy * dy + dz * dz);
  return 0;
}

int beacon_vector2str_array(const std::vector<ids::BeaconStation>& beacons,
    char**& strarray) noexcept
{
  strarray = new char*[beacons.size()];
  for (int i = 0; i < (int)beacons.size(); i++) {
    strarray[i] = new char[4];
#ifdef __GNUC__
    #ifndef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
#endif
#endif
    std::strncpy(strarray[i], beacons[i].m_station_id, 4);
#ifdef __GNUC__
    #ifndef __clang__
#pragma GCC diagnostic pop
#endif
#endif
  }
  return 0;
}

int beacon_shift_factor(const char* beacon_id, const ids::DorisObsRinex* rnx,
    int& k) noexcept
{
  auto it = std::find_if(rnx->stations().begin(), rnx->stations().end(),
      [&](const ids::BeaconStation& b) {
        return !std::strncmp(beacon_id, b.m_internal_code, 3);
      });

  if (it == rnx->stations().end())
    return 1;

  k = it->m_shift_factor;
  return 0;
}

auto match_beacon(const char* beacon_id,
    std::vector<LatestBeaconData>& prdata) noexcept
{
  auto it = std::find_if(prdata.begin(), prdata.end(),
      [&](const LatestBeaconData& db) {
        const char* str = db.obs.m_beacon_id;
        return !std::strncmp(beacon_id, str, 3);
      });
  return it;
}

double iono_correction(double l2ghz, double l400mhz) noexcept
{
  constexpr double gamma = ids::GAMMA_FACTOR;
  constexpr double sqrt_g = ids::GAMMA_FACTOR_SQRT;
  constexpr double denom = gamma - 1e0;
  return (l2ghz - sqrt_g * l400mhz) / denom;
}

// Starec[BC] PCO in mm
void starec_iono_free_pco(double* pco) noexcept
{
  using ids::AntennaOffset;
  using ids::GroundAntennaType;

  pco[0] = AntennaOffset<GroundAntennaType::Starec_B, 1>::dnorth();
  pco[1] = AntennaOffset<GroundAntennaType::Starec_B, 1>::deast();
  pco[2] = AntennaOffset<GroundAntennaType::Starec_B, 1>::dup();
  const double nl2 = AntennaOffset<GroundAntennaType::Starec_B, 2>::dnorth();
  const double el2 = AntennaOffset<GroundAntennaType::Starec_B, 2>::deast();
  const double ul2 = AntennaOffset<GroundAntennaType::Starec_B, 2>::dup();

  double dn1 = nl2 - pco[0];
  double de1 = el2 - pco[1];
  double du1 = ul2 - pco[2];

  constexpr double denom = ids::GAMMA_FACTOR - 1e0;
  pco[0] += dn1 / denom;
  pco[1] += de1 / denom;
  pco[2] += du1 / denom;
}

double iono_free_pcv(const ids::BeaconStation* beacon, double z, int& error) noexcept
{
  using ids::GroundAntennaType;

  GroundAntennaType type = beacon->type();

  double pcv_l1 = 0e0, pcv_l2 = 0e0;
  int err1 = 1, err2 = 1;
  switch (type) {
  case ids::GroundAntennaType::Starec_B:
    pcv_l1 = ids::AntennaOffset<GroundAntennaType::Starec_B, 1>::pcv(z, err1);
    pcv_l2 = ids::AntennaOffset<GroundAntennaType::Starec_B, 2>::pcv(z, err2);
    break;
  case ids::GroundAntennaType::Starec_C:
    pcv_l1 = ids::AntennaOffset<GroundAntennaType::Starec_C, 1>::pcv(z, err1);
    pcv_l2 = ids::AntennaOffset<GroundAntennaType::Starec_C, 2>::pcv(z, err2);
    break;
  case ids::GroundAntennaType::Alcatel:
    pcv_l1 = ids::AntennaOffset<GroundAntennaType::Alcatel, 1>::pcv(z, err1);
    pcv_l2 = ids::AntennaOffset<GroundAntennaType::Alcatel, 2>::pcv(z, err2);
  }
  error = err1 + err2;
  assert(!error);

  constexpr double denom = ids::GAMMA_FACTOR - 1e0;
  return (pcv_l2 - pcv_l1) / denom;
}

int main(int argc, char* argv[])
{
  if (argc != 5) {
    fprintf(stderr,
        "\nUsage: %s [RINEX_FILE] [SP3 FILE] [SINEX] [GTP3_GRID]\n",
        argv[0]);
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
  double xyz[3] = { 0 }, dxyz[3] = { 0 };

  // read the gpt3 5x5 grid file (so that we don't have to parse it every time)
  dso::gpt3::gpt3_grid gpt3_grid(dso::gpt3::gpt3_grid_attributes<dso::gpt3::Gpt3Grid::grid5x5>::num_lines);
  if (dso::gpt3::parse_gpt3_grid(argv[4], &gpt3_grid)) {
    fprintf(stderr, "[ERROR] Failed to parse the gpt3 grid file!\n");
    return 123;
  }

  // DORIS RINEX instance
  ids::DorisObsRinex rnx(argv[1]);
  rnx.print_metadata();

  // latest data records from all beacons
  std::vector<LatestBeaconData> prv_data;
  prv_data.reserve(rnx.stations().size());

  // index of the 2GHz phase measurement
  int l1_idx = rnx.get_observation_code_index(
      ids::ObservationCode { ids::ObservationType::phase, 1 });
  if (l1_idx < 0) {
    fprintf(stderr,
        "[ERROR] Failed to find ObservationType L1 in RINEX\'s observation "
        "types vector! (traceback: %s)\n",
        __func__);
    return 1;
  }

  // index of the 400MHz phase measurement (need for iono-free reduction)
  int l2_idx = rnx.get_observation_code_index(
      ids::ObservationCode { ids::ObservationType::phase, 2 });
  if (l2_idx < 0) {
    fprintf(stderr,
        "[ERROR] Failed to find ObservationType L2 in RINEX\'s observation "
        "types vector! (traceback: %s)\n",
        __func__);
    return 1;
  }

  // index of the F measurement (relative frequency offset)
  int f_idx = rnx.get_observation_code_index(
      ids::ObservationCode { ids::ObservationType::frequency_offset });
  if (f_idx < 0) {
    fprintf(stderr,
        "[ERROR] Failed to find ObservationType F in RINEX\'s observation "
        "types vector! (traceback: %s)\n",
        __func__);
    return 1;
  }

  // index of the P measurement (pressure)
  int p_idx = rnx.get_observation_code_index(
      ids::ObservationCode { ids::ObservationType::ground_pressure });
  if (p_idx < 0) {
    fprintf(stderr,
        "[ERROR] Failed to find ObservationType P in RINEX\'s observation "
        "types vector! (traceback: %s)\n",
        __func__);
    return 1;
  }

  // use the SINEX file to extract beacon coordinates (for reference time of
  // RINEX)
  ids::BeaconCoordinates* crd = new ids::BeaconCoordinates[rnx.stations().size()];
  char** sites_array = nullptr;
  beacon_vector2str_array(rnx.stations(), sites_array);
  if (ids::extrapolate_sinex_coordinates(argv[3], sites_array,
          rnx.stations().size(),
          rnx.ref_datetime(), crd)) {
    fprintf(stderr, "Failed to extrapolate beacon coordinates!\n");
    return 9;
  }
  // dealocate memory
  for (int i = 0; i < (int)rnx.stations().size(); i++)
    delete[] sites_array[i];
  delete[] sites_array;

  // transform coordinates to ellipsoidal; will be of use ...
  ids::BeaconCoordinates* crd_ell = new ids::BeaconCoordinates[rnx.stations().size()];
  for (int i = 0; i < (int)rnx.stations().size(); i++) {
    dso::car2ell<dso::ellipsoid::grs80>(
        crd[i].x, crd[i].y, crd[i].z, crd_ell[i].x, crd_ell[i].y, crd_ell[i].z);
  }

  // starec Iono-free PCO (do we have alcatele antennas?)
  double ifpco[3];
  starec_iono_free_pco(ifpco);
  for (int i = 0; i < (int)rnx.stations().size(); i++) {
    if (rnx.stations()[i].m_station_id[3] == 'A') {
      fprintf(stderr, "[ERROR] Non-STAREC antenna found ... need more code!\n");
      return 2;
    }
  }

  // just a counter for number of measurements per beacon ... to have a summary
  // at the end (and number of unconsidered obs)
  std::vector<std::pair<long, long>> beacon_obs(rnx.stations().size(), { 0L, 0L });

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
    sv_intrp.interpolate_at(tnow, xyz, dxyz);

    auto beaconobs_set = it.cblock.begin(); // Current beacons's observation set
    while (beaconobs_set != it.cblock.end()) {
      // internal 3char name of beacon
      const char* beacon_internal_id = beaconobs_set->m_beacon_id;

      // get this out of the way first; record that we have a new obs for this
      // beacon
      update_beacon_obs(rnx.beacon_internal_id2id(beacon_internal_id),
          beacon_obs, rnx.stations());

      // get shift factor and nominal frequency (Hz)
      int k;
      if (beacon_shift_factor(beacon_internal_id, it.rnx, k)) {
        fprintf(stderr,
            "[ERROR] Failed to match beacon shift factor for beacon %.3s\n",
            beacon_internal_id);
        return 101;
      }
      ids::beacon_nominal_frequency(k, nf2g, nf4m);

      // sv-beacon distance
      double rho;
      if (beacon_sv_distance(rnx.beacon_internal_id2id(beacon_internal_id), crd,
              rnx.stations().size(), xyz, rho))
        return 120;

      // get the iono correction
      double iono_cor = iono_correction(beaconobs_set->m_values[l1_idx].m_value,
          beaconobs_set->m_values[l2_idx].m_value);

      // check if the 2GHz observation is marked for discontinuity
      int discontinuity = (beaconobs_set->m_values[l1_idx].m_flag2 == '1') + (beaconobs_set->m_values[l2_idx].m_flag2 == '1');

      // get doppler count (if previous obs exists)
      auto beac_it = match_beacon(beacon_internal_id, prv_data);
      if (beac_it == prv_data.end()) { // no previous obs; just update the array
        if (!new_arc_on_discontinuity || (new_arc_on_discontinuity && !discontinuity))
          prv_data.emplace_back(
              LatestBeaconData { *beaconobs_set, tnow, rho, iono_cor });
        else
          update_missed_beacon_obs(rnx.beacon_internal_id2id(beacon_internal_id),
              beacon_obs, rnx.stations());

        // previous obs exist but this is a discontinuity; start new arc
      } else if (beac_it != prv_data.end() && ((new_arc_on_discontinuity && discontinuity))) {
        beac_it->start_new_arc = 1;
        update_missed_beacon_obs(rnx.beacon_internal_id2id(beacon_internal_id),
            beacon_obs, rnx.stations());

        // previous obs exist and is marked with new arc; re-ignite the arc, add
        // this as previous measurement
      } else if (beac_it != prv_data.end() && ((new_arc_on_discontinuity && !discontinuity && beac_it->start_new_arc))) {
        beac_it->start_new_arc = 0;
        beac_it->obs = *beaconobs_set;
        beac_it->t = tnow;
        // double rho0 = beac_it->sv_beacon_rho;
        beac_it->sv_beacon_rho = rho;
        beac_it->l2gh_iono_cor = iono_cor;

      } else {
        // previous data exists! get Ndop and update array (apply iono
        // correction if needed)
        double pvalue = beac_it->obs.m_values[l1_idx].m_value + apply_iono_correction * beac_it->l2gh_iono_cor;
        double cvalue = beaconobs_set->m_values[l1_idx].m_value + apply_iono_correction * iono_cor;
        double Ndop = cvalue - pvalue;

        // get delta time (receiver time)
        dso::nanoseconds Dt_nsec = tnow.delta_sec(beac_it->t);
        double Dtsec = Dt_nsec.to_fractional_seconds();

        // get the relative frequency offset (Df/f0)
        double Df0 = beaconobs_set->m_values[f_idx].m_value;

        // get the tropospheric correction ...
        ids::BeaconCoordinates* bpos_xyz = // ptr to the beacon's cartesian coordinates
            get_beacon_crd(rnx.beacon_internal_id2id(beacon_internal_id), crd,
                rnx.stations().size());
        assert(bpos_xyz);
        int crd_offset = bpos_xyz - crd;
        assert(crd_offset >= 0 && crd_offset < (int)rnx.stations().size());
        ids::BeaconCoordinates* bpos_ell = crd_ell + crd_offset; // ptr to the beacon's ellipsoidal coordinates

        // call gpt3 to get coefficients
        dso::gpt3_result g3out;
        if (dso::gpt3_fast(tnow, &bpos_ell->x, &bpos_ell->y, &bpos_ell->z, 1,
                0, &gpt3_grid, &g3out)) {
          fprintf(stderr, "[ERROR] Failed to compute gpt3!\n");
          return 20;
        }

        // compute mapping functions (need zenith distance first)
        double north, east, up;
        dso::car2top<dso::ellipsoid::grs80>(xyz[0], xyz[1], xyz[2], bpos_xyz->x,
            bpos_xyz->y, bpos_xyz->z, north,
            east, up);
        double distance, azimouth, zenith;
        dso::top2daz(north, east, up, distance, azimouth, zenith);
#ifdef dEBUG
        double elevation = std::atan2(up, std::sqrt(north * north + east * east));
        assert(zenith >= 0e0 && zenith <= M_PI / 2e0 && std::abs(dso::rad2deg(zenith - (M_PI / 2 - elevation) < 0.1e0)));
#endif
        double mfh, mfw;
        dso::vmf3(g3out.ah, g3out.aw, tnow, bpos_ell->x, bpos_ell->y, zenith,
            mfh, mfw);

        // ZTD hydrostatic
        double gpres_hpa = beaconobs_set->m_values[p_idx].m_value; // pressure value from RINEX
        double ztd_hydro = dso::saasthyd(gpres_hpa, bpos_ell->x,
            bpos_ell->z); // ZTD hydrostatic

        // TZD wet
        double ztd_wet = dso::asknewet(g3out.e, g3out.Tm, g3out.la);

        // get PCV (in mm)
        int err1 = 0, err2 = 0;
        double pcv_l1 = ids::AntennaOffset<ids::GroundAntennaType::Starec_B, 1>::pcv(zenith,
            err1);
        double pcv_l2 = ids::AntennaOffset<ids::GroundAntennaType::Starec_B, 2>::pcv(zenith,
            err2);
        assert(!(err1 + err2));

        // update the data array
        beac_it->obs = *beaconobs_set;
        beac_it->t = tnow;
        double rho0 = beac_it->sv_beacon_rho;
        beac_it->sv_beacon_rho = rho;
        beac_it->l2gh_iono_cor = iono_cor;

        printf("%.3s %.5f %.5f %.5f %.8f %.5f %.8f t: %.2f %.6f %.6f %.6f %.6f "
               ":t %.5f %.5f %.6f %.6f [%c%c]\n",
            beacon_internal_id, tnow.as_mjd(), nf2g, Ndop, Dtsec, Df0,
            (rho - rho0) / Dtsec, dso::rad2deg(zenith), ztd_hydro, mfh,
            ztd_wet, mfw, ztd_hydro * mfh + ztd_wet * mfw, iono_cor,
            pcv_l1, pcv_l2,
            beaconobs_set->m_values[l1_idx].m_flag1,
            beaconobs_set->m_values[l1_idx].m_flag2);
      }
      ++beaconobs_set;
    }
  }

  unsigned long all_obs = 0, flagged_obs = 0;
  printf("Read lines; ended with %d\n", error);
  printf("Number of epochs per beacon:\n");
  int index = 0;
  for (const auto& b : rnx.stations()) {
    double precent_obs_removed = (double)beacon_obs[index].second * 100 / (double)beacon_obs[index].first;
    printf(
        "%.4s (%.3s) %6ld/%6ld missed: %.1f%%\n", b.m_station_id, b.m_internal_code,
        beacon_obs[index].first, beacon_obs[index].second, precent_obs_removed);
    all_obs += beacon_obs[index].first;
    flagged_obs += beacon_obs[index].second;
    ++index;
  }
  printf("Total number of L2GHz obs=%lu, flagged=%lu, percent flagged=%.1f%%\n",
      all_obs, flagged_obs, flagged_obs * 100 / (double)all_obs);

  printf("Iono-free PCO: n=%.3f e=%.3f u=%.3f\n", ifpco[0], ifpco[1], ifpco[2]);

  return 0;
}
