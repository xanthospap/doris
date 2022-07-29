#ifndef __DSO_DORIS_VAR_UTILS_HPP__
#define __DSO_DORIS_VAR_UTILS_HPP__

#include "datetime/dtcalendar.hpp"
#include "var/celestrak.hpp"

namespace dso {

int get_CelesTrack_flux_data(
    const dso::modified_julian_day mjd, const char *fncsv,
    dso::utils::celestrak::details::CelestTrakSWFlux &flux_data) noexcept {
  return dso::utils::celestrak::details::parse_csv_for_date(mjd, fncsv,
                                                             flux_data);
}

template <typename T>
int get_CelesTrack_flux_data(
    const dso::datetime<T> &t, const char *fncsv,
    dso::utils::celestrak::details::CelestTrakSWFlux &flux_data) noexcept {
  return dso::utils::celestrak::details::parse_csv_for_date(t.as_mjd(), fncsv,
                                                             flux_data);
}

} // namespace dso

#endif
