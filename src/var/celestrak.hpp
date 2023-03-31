#ifndef __DSO__CELESTRACK_VAR_DATA_UTILS_HPP__
#define __DSO__CELESTRACK_VAR_DATA_UTILS_HPP__

#include "datetime/dtcalendar.hpp"
#include <fstream> // for pos_type

namespace dso::utils::celestrak::details {

constexpr const double MissingSwData = std::numeric_limits<double>::min();

/// See https://celestrak.org/SpaceData/
/// Flux related data from Celestrack Space Wether (CSV) data
/// see https://celestrak.org/SpaceData/SpaceWx-format.php
struct CelestTrakSWFlux {
  dso::modified_julian_day mjd_;
  double f107Obs{MissingSwData};        ///< F10.7_OBS -> 25
  double f107Adj{MissingSwData};           ///< F10.7_ADJ -> 26
  double f107ObsC81{MissingSwData};        ///< F10.7_OBS_CENTER81 -> 28
  double f107AdjC81{MissingSwData};        ///< F10.7_ADJ_CENTER81 -> 30
  double f107ObsL81{MissingSwData};        ///< F10.7_OBS_LAST81 -> 29
  double f107AdjL81{MissingSwData};        ///< F10.7_ADJ_LAST81 -> 31
  double ApDailyAverage{MissingSwData}; ///< Arithmetic average of the 8 Ap
                                        ///< indices for the day. -> 21
  int ApIndexes[8];                     ///< 3-hour Ap indexes -> [13-20]
  int KpIndexes[8];                     ///< 3-hour Kp indexes -> [3-10]
  double KpSum{MissingSwData};          ///< Sum of the 8 Kp indices for the day
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

// TODO the following two functions are obsolete and should be deleted
/// @brief Parse a SW CSV file and return Flux info for a given number of days,
///        centered around given date.
/// The function will parse a (what is expected to be) SW CSV file and return
/// individual CelestTrakSWFlux instances, for
/// * mjd - days_before
/// * mjd - (days_before - 1)
/// * mjd - (days_before - ...)
/// * mjd
/// * mjd + 1
/// * mjd + ...
/// * mjd + days_after
/// SW CSV file(s), can be downloaded from CelesTrak, see
/// https://celestrak.org/SpaceData/
/// @param[in] mjd Center MJD for requested Flux data
/// @param[in] fncsv The filename of the SW SCV data file
/// @param[out] flux_data An array of CelestTrakSWFlux instances;
///             at output they hold respective data for the dates:
///             [mjd-days_before, ... , mjd-1, mjd, mjd+1, mjd+days_after]
///             It must be of size (at least) : days_before + 1 + days_after
/// @param[out] fpos Position of tellg() on the input file, after retrieving
///             the last line. This can be later used to easily get the next
///             record, for the next day
/// @param[in] days_before Number of days to be collected, before mjd
/// @param[in] days_after Number of days to be collected, after mjd
/// @return Anything other than 0 denotes an error. In this case, the
///         resulting flux_data may be erronuous!
int parse_csv_for_date(dso::modified_julian_day mjd, const char *fncsv,
                       CelestTrakSWFlux *flux_data,
                       std::ifstream::pos_type &fpos, int days_before = 4,
                       int days_after = 0) noexcept;

int get_next(dso::modified_julian_day mjd, const char *fncsv,
             CelestTrakSWFlux *flux_data, int flux_data_sz,
             std::ifstream::pos_type &fpos) noexcept;

} // namespace dso::utils::celestrak::details
#endif
