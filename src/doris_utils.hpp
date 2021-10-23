#ifndef __DORIS_SYS_UTILS_HPP__
#define __DORIS_SYS_UTILS_HPP__

#include "doris_system_info.hpp"
#include "ggdatetime/dtcalendar.hpp"
#include <type_traits>

namespace ids {

/// @brief Extrapolate site coordinates using information from a SINEX file
/// This function is mainly designed to work with DPOD SINEX files (see
/// https://ids-doris.org/analysis-coordination/combination/dpod.html), though
/// it could work with any SINEX file that has SOLUTION/ESTIMATE block holding
/// STA[XYZ] and VEL[XYZ] parameter entries.
/// @param[in] snx_in The filename of the SINEX file to use; this file should
///            have a SOLUTION/ESTIMATE block; we will seek the STA[XYZ] and
///            VEL[XYZ] parameters for the requested sites
/// @param[in] station_ids Array of strings, aka site id (4-char) strings; for
///            each one we will try to extrapolate its coordinates to the
///            epoch requested (e.g {..., "DIOA", ...})
/// @param[in] num_sites Number of sites in the station_ids array
/// @param[in] t The epoch we want to extrapolate coordinates at
/// @return An integer; Anything other than 0 denotes an error
int extrapolate_sinex_coordinates(
    const char *snx_fn, char **site_ids, int num_sites,
    const dso::datetime<dso::microseconds> &t) noexcept;

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
    const char *snx_fn, char **site_ids, int num_sites,
    const dso::datetime<T> &t) noexcept {
  dso::datetime<dso::microseconds> t_micro =
      t.template cast_to<dso::microseconds>();
  return extrapolate_sinex_coordinates(snx_fn, site_ids, num_sites, t_micro);
}
} // namespace ids

#endif