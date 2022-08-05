#ifndef __DORIS_SYS_UTILS_HPP__
#define __DORIS_SYS_UTILS_HPP__

#include "datetime/dtcalendar.hpp"
#include "doris_system_info.hpp"
#include <vector>
#include <type_traits>

namespace dso {

/// @brief This structure is meant to hold beacon coordinates and std.
/// deviations
struct BeaconCoordinates {
  char id[4];
  ///< x,y,z (cartesian) coordinate components in [m]
  double x, y, z;
  ///< std. deviation values of the above components in [m^2]
  double xs, ys, zs;
}; // BeaconCoordinates

/// @brief Extrapolate site coordinates using information from a SINEX file
/// This function is mainly designed to work with DPOD SINEX files (see
/// https://ids-doris.org/analysis-coordination/combination/dpod.html), though
/// it could work with any SINEX file that has SOLUTION/ESTIMATE block holding
/// STA[XYZ] and VEL[XYZ] parameter entries.
///
/// @warning Users should use the wrapper function that works with std::vector
/// (instead of pointers), see below.
///
/// @param[in] snx_in The filename of the SINEX file to use; this file should
///            have a SOLUTION/ESTIMATE block; we will seek the STA[XYZ] and
///            VEL[XYZ] parameters for the requested sites
/// @param[in] station_ids Array of strings, aka site id (4-char) strings; for
///            each one we will try to extrapolate its coordinates to the
///            epoch requested (e.g {..., "DIOA", ...})
/// @param[in] num_sites Number of sites in the station_ids array
/// @param[in] t The epoch we want to extrapolate coordinates at
/// @param[out] BeaconCoordinates An array of BeaconCoordinates to hold
///            results, aka station id's and extrapolated coordinates at t.
///            Its size must be at least num_sites. If it is larger, or some
///            of the requested sites were not found in the SINEX file, the
///            (any) remaining elements of BeaconCoordinates amy contain
///            garbage! To get the actual number of sites extrapolated, use
///            sites_extrapolated; all remaining elements of BeaconCoordinates
///            are garbage.
/// @param[out] sites_extrapolated Number of sites for which the coordinates
///            were extrapolated and written to BeaconCoordinates
/// @param[in] missing_site_is_error If set to true, the function will
///            stop with error if any of the sites were not found in the SINEX
///            file.
/// @return An integer; Anything other than 0 denotes an error
int extrapolate_sinex_coordinates(const char *snx_fn, char **site_ids,
                                  int num_sites,
                                  const dso::datetime<dso::microseconds> &t,
                                  BeaconCoordinates *results,
                                  int &sites_extrapolated,
                                  bool missing_site_is_error = true) noexcept;

/// @brief Extrapolate site coordinates using information from a SINEX file
/// This function is mainly designed to work with DPOD SINEX files (see
/// https://ids-doris.org/analysis-coordination/combination/dpod.html), though
/// it could work with any SINEX file that has SOLUTION/ESTIMATE block holding
/// STA[XYZ] and VEL[XYZ] parameter entries.
///
/// @param[in] snx_in The filename of the SINEX file to use; this file should
///            have a SOLUTION/ESTIMATE block; we will seek the STA[XYZ] and
///            VEL[XYZ] parameters for the requested sites
/// @param[in] station_ids Array of dso::BeaconStation which we will try to
///            extrapolate coordinates to the epoch requested. To identify
///            the dso::BeaconStation to a site written in the SINEX file, we
///            will use the dso::BeaconStation::m_station_id member variable,
///            aka it's 4-char ID.
/// @param[in] t The epoch we want to extrapolate coordinates at
/// @param[out] BeaconCoordinates An array of BeaconCoordinates to hold
///            results, aka station id's and extrapolated coordinates at t.
/// @param[in] missing_site_is_error If set to true, the function will
///            stop with error if any of the sites were not found in the SINEX
///            file.
/// @return An integer; Anything other than 0 denotes an error
int extrapolate_sinex_coordinates(
    const char *snxfn, const std::vector<dso::BeaconStation> &station_ids,
    const dso::datetime<dso::microseconds> &t,
    std::vector<dso::BeaconCoordinates> &result_array,
    bool missing_site_is_error = true) noexcept;

/// @brief Overload of the above, for any datetime<T> type.
#if __cplusplus >= 202002L
template <typename T>
requires(T::is_of_sec_type && !std::is_same_v<T, dso::microseconds>)
#else
template <class T,
          typename = std::enable_if_t<(T::is_of_sec_type) &&
                                      (!std::is_same_v<T, dso::microseconds>)>>
#endif
    int extrapolate_sinex_coordinates(
        const char *snxfn, const std::vector<dso::BeaconStation> &station_ids,
        const dso::datetime<T> &t,
        std::vector<dso::BeaconCoordinates> &result_array,
        bool missing_site_is_error) noexcept {
  dso::datetime<dso::microseconds> t_micro =
      t.template cast_to<dso::microseconds>();
  return extrapolate_sinex_coordinates(snxfn, station_ids, t_micro,
                                       result_array, missing_site_is_error);
}
} // namespace ids

#endif
