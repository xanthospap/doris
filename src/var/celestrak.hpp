#ifndef __DSO__CELESTRACK_VAR_DATA_UTILS_HPP__
#define __DSO__CELESTRACK_VAR_DATA_UTILS_HPP__

#include "datetime/dtcalendar.hpp"

namespace dso::utils::celestrak::details {

constexpr const double MissingSwData = std::numeric_limits<double>::min();

// See
// https://celestrak.org/SpaceData/
/// Flux related data from Celestrack Space Wether (CSV) data
/// see https://celestrak.org/SpaceData/SpaceWx-format.php
struct CelestTrakSWFlux {
  dso::modified_julian_day mjd_;
  double f107Obs{MissingSwData},        ///< F10.7_OBS -> 25
      f107Adj{MissingSwData},           ///< F10.7_ADJ -> 26
      f107ObsC81{MissingSwData},        ///< F10.7_OBS_CENTER81 -> 28
      f107ObsL81{MissingSwData},        ///< F10.7_OBS_LAST81 -> 29
      f107AdjC81{MissingSwData};        ///< F10.7_ADJ_CENTER81 -> 30
  double ApDailyAverage{MissingSwData}; ///< Arithmetic average of the 8 Ap
                                        ///< indices for the day. -> 21
  double ApIndexes[8];                  ///< 3-hour Ap indexes -> [13-20]
  char flag[5] = {'\0'};                ///< Flux Qualifier
};

/// @brief Resolve a SW CSV line record to a CelestTrakSWFlux instance
/// @param[in] line The SW CSV line record to resolve (see class documentation
///                 for data format)
/// @param[out] flux_data A CelestTrakSWFlux info holding parsed data from the
///                 input line.
/// @return Anything other than 0 denotes an error. In this case, do not use
///         the flux_data.
int resolve_csv_line_records(const char *line,
                             CelestTrakSWFlux &flux_data) noexcept;

/// @brief Resolve a SW CSV line record's date. These line only hold date
///        info, no time at all.
/// @param[in] line The SW CSV line record to resolve (see class documentation
///                 for data format)
/// @param[out] mjd The record date as MJD
/// @return Anything other than 0 denotes an error. In this case, the resulting
///         mjd may be erronuous!
int resolve_csv_line_date(const char *line,
                          dso::modified_julian_day &mjd) noexcept;

/// @brief Parse a SW CSV file and return Flux info for three days, the middle
///        one being mjd.
/// The function will parse a (what is expected to be) SW CSV file and return
/// three individual CelestTrakSWFlux instances, for
/// * mjd - 1
/// * mjd, and
/// * mjd + 1
/// aka three days, centered at given mjd.
/// SW CSV file(s), can be downloaded from CelesTrak, see
/// https://celestrak.org/SpaceData/
/// @param[in] mjd Center MJD for requested Flux data
/// @param[in] fncsv The filename of the SW SCV data file
/// @param[out] flux_data An array of (at least) 3 CelestTrakSWFlux instances;
///             at output they hold respective data for the dates:
///             [mjd-1, mjd, mjd+1]
/// @return Anything other than 0 denotes an error. In this case, the resulting
///         flux_data may be erronuous!
int parse_csv_for_date(dso::modified_julian_day mjd, const char *fncsv,
                       CelestTrakSWFlux *flux_data) noexcept;
} // namespace dso::utils::celestrak::details
#endif
