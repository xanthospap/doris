#ifndef __DORIS_SYS_UTILS_HPP__
#define __DORIS_SYS_UTILS_HPP__

#include "doris_system_info.hpp"
#include "ggdatetime/dtcalendar.hpp"

namespace ids {
#if __cplusplus >= 202002L
template <gconcepts::is_sec_dt T>
#else
template <class T, typename = std::enable_if_t<
    std::enable_if_t<T::is_of_sec_type> && 
    std::is_same_t<T, dso::microseconds>>>
    >
#endif
int extrapolate_sinex_coordinates(const char *snx_fn, char **station_ids,
                                  int num_stations,
                                  const dso::datetime<T> &t) noexcept {
  dso::datetime<dso::microseconds> t_micro = t. template cast_to<dso::microseconds>();
  printf(">> casting datetime to microsecnds\n");
  return extrapolate_sinex_coordinates(snx_fn, station_ids, num_stations,
                                       t_micro);
}

int extrapolate_sinex_coordinates(const char *snx_fn, char **station_ids,
                                  int num_stations,
                                  const dso::datetime<dso::microseconds> &t) noexcept;
} // namespace ids

#endif