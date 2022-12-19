#include "doris_system_info.hpp"
#include "doris_utils.hpp"
#include "sinex/sinex.hpp"
#include <cstring>

int extrapolate_coordinates(
    const char *siteid, std::vector<dso::sinex::SolutionEstimate> &estimates,
    const dso::datetime<dso::microseconds> &t, double *pos) noexcept {

  const char *components[] = {"STAX", "STAY", "STAZ"};
  for (int i = 0; i < 3; i++) {
    const char *sta = components[i];
    char vel[4] = "VEL";
    vel[3] = sta[3];
    auto pos_sx =
        std::find_if(estimates.begin(), estimates.end(),
                     [&](const dso::sinex::SolutionEstimate &est) {
                       return !std::strncmp(est.m_site_code, siteid, 4) &&
                              !std::strncmp(est.m_param_type, sta, 4);
                     });
    auto pos_vx =
        std::find_if(estimates.begin(), estimates.end(),
                     [&](const dso::sinex::SolutionEstimate &est) {
                       return !std::strncmp(est.m_site_code, siteid, 4) &&
                              !std::strncmp(est.m_param_type, vel, 4);
                     });

    // check that we have position and velocity
    if (pos_sx == estimates.end() || pos_vx == estimates.end()) {
      fprintf(stderr,
              "[ERROR] %s/%s missing from collected estimates (site %s); "
              "cannot extrapolate coordinates (*traceback: %s)\n",
              sta, vel, siteid, __func__);
      return 1;
    }

    // check that the units are ok
    if (std::strncmp(pos_sx->m_units, "m ", 2) ||
        std::strncmp(pos_vx->m_units, "m/y", 3)) {
      fprintf(stderr,
              "[ERROR] Invalid units for %s/%s estimates (site %s); cannot "
              "extrapolate coordinates (*traceback: %s)\n",
              sta, vel, siteid, __func__);
      return 1;
    }

    // delta time in years
    double mjd_ref = pos_sx->m_epoch.as_mjd();
    double mjd_at = t.as_mjd();
    double dt = (mjd_at - mjd_ref) / 365.25e0;

    // extrapolate the coordinate
    pos[i] = pos_sx->m_estimate + pos_vx->m_estimate * dt;
  }
  return 0;
}

int dso::extrapolate_sinex_coordinates(
    const char *snxfn, const std::vector<dso::BeaconStation> &beacons,
    const dso::datetime<dso::microseconds> &t,
    std::vector<dso::BeaconCoordinates> &result_array,
    bool missing_site_is_error) noexcept {
  // allocate and fill the list of site id's
  int num_sites = beacons.size();
  char **sites = new char *[num_sites];
  for (int i = 0; i < num_sites; i++) {
    sites[i] = new char[5];
    std::memset(sites[i], 0, 5);
    std::strncpy(sites[i], beacons[i].m_station_id, 4);
  }

  // check that the result_array has enough memmory to hold results
  // also, clear it up
  result_array.clear();
  // We are going to use the raw memory of the std::vector, so do not remove
  // the following line; if needed, resize at the end
  result_array.resize(beacons.size());

  // call the core function to extrapolate the coordinates
  int sites_found;
  int error = dso::extrapolate_sinex_coordinates(
      snxfn, sites, num_sites, t, result_array.data(), sites_found,
      missing_site_is_error);

  // successeful or not, we do not need the sites array anymore
  // deallocate memory
  for (int i = 0; i < num_sites; i++)
    delete[] sites[i];
  delete[] sites;

  // resize the output vector so we don't use garbage values
  result_array.resize(sites_found);

  // all done
  return error;
}

int dso::extrapolate_sinex_coordinates(
    const char *snx_fn, char **station_ids, int num_stations,
    const dso::datetime<dso::microseconds> &t,
    dso::BeaconCoordinates *result_array, int &result_size,
    bool missing_site_is_error) noexcept {

  // clear result size
  result_size = 0;

  // a vector to hold indexes of station_ids that are missing from the SINEX
  // file
  std::vector<int> missing_indexes;

  // initialize the sinex instance
  dso::Sinex snx(snx_fn);

  // a vector of SiteId to hold results
  std::vector<dso::sinex::SiteId> site_vec;

  // parse the block SITE/ID to collect info for the given sites
  int error = snx.parse_block_site_id(site_vec, num_stations, station_ids);
  if (error || (int)site_vec.size() != num_stations) {
    if (error) {
      fprintf(
          stderr,
          "[ERROR] Failed to collect info for given sites; SINEX file is %s "
          "(traceback: %s)\n",
          snx.filename().c_str(), __func__);
      return 1;
    } else {
      // report the missing station(s)
      for (int i = 0; i < num_stations; i++) {
        const char *sid = station_ids[i];
        const auto it =
            std::find_if(site_vec.begin(), site_vec.end(),
                         [&](const dso::sinex::SiteId &rid) {
                           return !std::strncmp(sid, rid.m_site_code, 4);
                         });
        if (it == site_vec.end()) {
          fprintf(stderr,
                  "[WARNING] Failed to find station %.4s in SINEX file %s\n",
                  sid, snx_fn);
          missing_indexes.push_back(i);
        }
      }
      if (missing_site_is_error)
        return 2;
    }
  }

  // first declare a vector of SolutionEstimate to hold results
  std::vector<dso::sinex::SolutionEstimate> est_vec;

  // get the solution estimates for the sites
  error = snx.parse_block_solution_estimate(est_vec, site_vec);
  if (error || (int)est_vec.size() < num_stations * 6) {
    fprintf(stderr,
            "[ERROR] Failed to collect estimates for given sites; SINEX file "
            "is %s "
            "(traceback: %s)\n",
            snx.filename().c_str(), __func__);
    return 1;
  }

  // for each of the stations, extrapolate the coordinate estimates per
  // component
  double pos[3];
  for (int i = 0; i < num_stations; i++) {
    const char *sid = station_ids[i];

    // only proceed, if station is not missing!
    if (std::vector<int>::iterator it;
        (it = std::find_if(missing_indexes.begin(), missing_indexes.end(),
                           [&](int idx) { return i == idx; })) ==
        missing_indexes.end()) {

      if (extrapolate_coordinates(sid, est_vec, t, pos)) {
        fprintf(stderr,
                "[ERROR] Failed to extrapolate coordinates for site %s from "
                "SINEX %s (traceback: %s)\n",
                sid, snx.filename().c_str(), __func__);
      } else {
        std::memcpy(result_array[result_size].id, sid, 4 * sizeof(char));
        result_array[result_size].x = pos[0];
        result_array[result_size].y = pos[1];
        result_array[result_size].z = pos[2];
        ++result_size;
      }
    }
  }

  return 0;
}
