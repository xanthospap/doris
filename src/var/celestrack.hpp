#ifndef __DSO__CELESTRACK_VAR_DATA_UTILS_HPP__
#define __DSO__CELESTRACK_VAR_DATA_UTILS_HPP__

#include "datetime/dtcalendar.hpp"

namespace dso::utils::celestrack::details {

constexpr const double MissingSwData = std::numeric_limits<double>::min();

// See
// https://celestrak.org/SpaceData/
/// Flux related data from Celestrack Space Wether (CSV) data
struct CelestTrackSWFlux {
  dso::modified_julian_day mjd_;
  double f107Obs{MissingSwData}, ///< F10.7_OBS -> 25
      f107Adj{MissingSwData},    ///< F10.7_ADJ -> 26
      f107ObsC81{MissingSwData}, ///< F10.7_OBS_CENTER81 -> 28
      f107ObsL81{MissingSwData}, ///< F10.7_OBS_LAST81 -> 29
      f107AdjC81{MissingSwData}; ///< F10.7_ADJ_CENTER81 -> 30
};

int resolve_csv_line_records(const char *line, double *data) noexcept;

int resolve_csv_line_date(const char *line,
                          dso::modified_julian_day &mjd) noexcept;

int parse_csv_for_date(dso::modified_julian_day mjd, const char *fncsv,
                       CelestTrackSWFlux &flux_data) noexcept;
} // namespace dso::utils::celestrack::details
#endif
