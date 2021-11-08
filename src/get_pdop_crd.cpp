#include "doris_utils.hpp"
#include "sinex.hpp"
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

int ids::extrapolate_sinex_coordinates(
    const char *snx_fn, char **station_ids, int num_stations,
    const dso::datetime<dso::microseconds> &t, ids::BeaconCoordinates *result_array) noexcept {
  
  // initialize the sinex instance
  dso::Sinex snx(snx_fn);

  // a vector of SiteId to hold results
  std::vector<dso::sinex::SiteId> site_vec;

  // parse the block SITE/ID to collect info for the given sites
  int error = snx.parse_block_site_id(site_vec, num_stations, station_ids);
  if (error || (int)site_vec.size() != num_stations) {
    fprintf(stderr,
            "[ERROR] Failed to collect info for given sites; SINEX file is %s "
            "(traceback: %s)\n",
            snx.filename().c_str(), __func__);
    return 1;
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
    if (extrapolate_coordinates(sid, est_vec, t, pos)) {
      fprintf(stderr,
              "[ERROR] Failed to extrapolate coordinates for site %s from "
              "SINEX %s (traceback: %s)\n",
              sid, snx.filename().c_str(), __func__);
      return 1;
    }
    std::memcpy(result_array[i].id, sid, 4);
    result_array[i].x = pos[0];
    result_array[i].y = pos[1];
    result_array[i].z = pos[2];
    // printf("%s %20.5f %20.5f %20.5f\n", sid, pos[0], pos[1], pos[2]);
  }

  return 0;
}