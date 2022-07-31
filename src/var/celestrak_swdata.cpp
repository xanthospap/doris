#include "celestrak.hpp"
#include <charconv>
#include <cstdio>
#include <cstring>
#include <datetime/dtfund.hpp>
#include <fstream>
#include <limits>
#include <system_error>

constexpr const int fn_error = -1;
constexpr const int parse_error = -2;
constexpr const int target_prior = 1;
constexpr const int eof_reached = 9;
constexpr const int MAX_SW_CHARS = 512;
constexpr const char header[] =
    "DATE,BSRN,ND,KP1,KP2,KP3,KP4,KP5,KP6,KP7,KP8,KP_SUM,AP1,AP2,AP3,AP4,AP5,"
    "AP6,AP7,AP8,AP_AVG,CP,C9,ISN,F10.7_OBS,F10.7_ADJ,F10.7_DATA_TYPE,F10.7_"
    "OBS_CENTER81,F10.7_OBS_LAST81,F10.7_ADJ_CENTER81,F10.7_ADJ_LAST81";

// F10.7_OBS -> 25
// F10.7_ADJ -> 26
// F10.7_OBS_CENTER81 -> 28
// F10.7_OBS_LAST81 -> 29
// F10.7_ADJ_CENTER81 -> 30

int dso::utils::celestrak::details::resolve_csv_line_records(
    const char *line,
    dso::utils::celestrak::details::CelestTrakSWFlux &flux_data) noexcept {

  // parse date
  if (dso::utils::celestrak::details::resolve_csv_line_date(line,
                                                            flux_data.mjd_)) {
    return 1;
  }

  const auto sz = std::strlen(line);
  const char *c = line;

  // find and skip 12 ',' chars
  int commas = 0;
  while (*c++ && commas < 12) {
    if (*c == ',') {
      ++commas;
    }
  }
  if (!*c) {
    return 2;
  }

  // last character in the string, the limit for conversions
  const char *last = line + sz - 1;

  // resolve the 8 hourly AP's (these are ints)
  int i = 0;
  while (i < 8) {
    auto pec = std::from_chars(c, last, flux_data.ApIndexes[i]);
    c = pec.ptr;
    if (*c++ != ',' || pec.ec != std::errc{}) {
      return 3 + i;
    }
    ++i;
  }

  // now ptr should point to the first char (after comma) of the next column
  // next column is daily average Ap
  auto pec = std::from_chars(c, last, flux_data.ApDailyAverage);
  if (pec.ec != std::errc{} || *pec.ptr++ != ',') {
    return 11;
  }
  c = pec.ptr;

  // we are now at column CP; skip and, and do so also for columns:
  // C9, ISN
  // Go three ',' forward
  commas = 0;
  while (*c++ && commas < 3) {
    if (*c == ',') {
      ++commas;
    }
  }
  if (!*c) {
    return 12;
  }

  // we are not at the F10.7_OBS column (25)
  pec = std::from_chars(c, last, flux_data.f107Obs);
  c = pec.ptr;
  if (pec.ec != std::errc{} || *c++ != ',') {
    return 13;
  }

  // we are not at the F10.7_ADJ
  pec = std::from_chars(c, last, flux_data.f107Adj);
  c = pec.ptr;
  if (pec.ec != std::errc{} || *c++ != ',') {
    return 14;
  }

  // we are now at the F10.7_DATA_TYPE; just copy three chars
  std::memcpy(flux_data.flag, c, 3 * sizeof(char));
  c += 3;
  if (*c++ != ',')
    return 15;

  // column F10.7_OBS_CENTER81
  pec = std::from_chars(++c, last, flux_data.f107ObsC81);
  c = pec.ptr;
  if (pec.ec != std::errc{} || *c++ != ',')
    return 16;

  // column F10.7_OBS_LAST81
  pec = std::from_chars(++c, last, flux_data.f107ObsL81);
  c = pec.ptr;
  if (pec.ec != std::errc{} || *c++ != ',')
    return 17;

  // column F10.7_ADJ_CENTER81
  pec = std::from_chars(++c, last, flux_data.f107AdjC81);
  c = pec.ptr;
  if (pec.ec != std::errc{} || *c++ != ',')
    return 18;

  // column F10.7_ADJ_LAST81
  // pec = std::from_chars(++c, last, flux_data.f107Adj);
  // c = pec.ptr;
  // if (pec.ec != std::errc{} || *c++ != ',')
  //  return 1;

  return 0;
}

int dso::utils::celestrak::details::resolve_csv_line_date(
    const char *line, dso::modified_julian_day &mjd) noexcept {
  // date format: 2017-01-10
  int year, month, day;
  // year ...
  auto pec = std::from_chars(line, line + 4, year);
  if (pec.ec != std::errc{})
    return 1;
  // month ...
  pec = std::from_chars(line + 5, line + 7, month);
  if (pec.ec != std::errc{})
    return 1;
  // day ...
  pec = std::from_chars(line + 8, line + 10, day);
  if (pec.ec != std::errc{})
    return 1;

  // to datetime ...
  const dso::datetime<dso::seconds> d(dso::year(year), dso::month(month),
                                      dso::day_of_month(day), dso::seconds(0));
  mjd = d.mjd();

  // all done
  return 0;
}

int dso::utils::celestrak::details::parse_csv_for_date(
    dso::modified_julian_day mjd, const char *fncsv,
    dso::utils::celestrak::details::CelestTrakSWFlux *flux_data,
    int days_before, int days_after) noexcept {

  // open the csv file
  std::ifstream fin(fncsv);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR@%s: Failed opening SW csv file %s\n", __func__,
            fncsv);
    return fn_error;
  }

  // set mjd's
  dso::modified_julian_day start_mjd(mjd.as_underlying_type() - days_before);
  dso::modified_julian_day end_mjd(mjd.as_underlying_type() + days_after);

  // read first line; should match the expected header
  char line[MAX_SW_CHARS];
  const auto hsz = std::strlen(header);
  if (!fin.getline(line, MAX_SW_CHARS)) {
    return fn_error;
    // may include whitespace chars ... no std::strcmp
  } else if (std::strncmp(line, header, hsz)) {
    fprintf(stderr, "ERROR@%s: Failed to match header line in SW csv file %s\n",
            __func__, fncsv);
    return fn_error;
  }

  dso::modified_julian_day cmjd;

  // read first data line; if record_date > start_mjd data, the target date is
  // not included in the file ...
  fin.getline(line, MAX_SW_CHARS);
  if (resolve_csv_line_date(line, cmjd)) {
    fprintf(stderr, "ERROR@%s: Failed parsing SW csv line %s\n", __func__,
            line);
    return parse_error;
  }
  if (cmjd > start_mjd)
    return target_prior;

  // keep on parsing/reading lines, untill we reach start_mjd or EOF
  int error = 0;
  while (cmjd != start_mjd && !error) {
    if (fin.getline(line, MAX_SW_CHARS))
      error = resolve_csv_line_date(line, cmjd);
    else
      error = 999;
  }

  // check for error/EOF
  if (error) {
    if (fin.eof())
      return eof_reached;
    return error;
  }

  int days_to_resolve = days_before + 1 + days_after;
  int cdi = 0;
  while (cdi < days_to_resolve) {
    // resolve current date data, store in flux_data[cdi]
    if ((error = resolve_csv_line_records(line, flux_data[cdi]))) {
      fprintf(stderr,
              "ERROR@%s: Failed resolving flux data from line %s (error=%d)\n",
              __func__, line, error);
      return parse_error;
    }
    // read in next line (should be next day)
    if (fin.getline(line, MAX_SW_CHARS)) {
      error = resolve_csv_line_date(line, cmjd);
    } else {
      return error;
    }
    if (cmjd != start_mjd + dso::modified_julian_day{cdi + 1}) {
      fprintf(stderr, "ERROR@%s Error parsing SW CVS file %s\n", __func__,
              fncsv);
      return parse_error;
    }
    ++cdi;
  }

  return 0;
}
