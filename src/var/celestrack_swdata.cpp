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

int dso::utils::celestrak::details::resolve_csv_line_records(const char *line,
                                                     double *data) noexcept {
  const auto sz = std::strlen(line);
  const char *c = line;

  // find and skip 24 ',' chars
  int commas = 0;
  while (*c++ && commas < 24) {
    if (*c == ',') {
      ++commas;
    }
  }
  if (!*c) {
    return 1;
  }

  // last character in the string, the limit for conversions
  const char *last = line + sz - 1;

  // resolve ...
  auto pec = std::from_chars(c, last, data[0]);
  if (pec.ec != std::errc{})
    return 1;

  // now ptr should point to the first, non-resolved char, aka next comma!
  c = pec.ptr;
  if (*c++ != ',') {
    return 1;
  }
  pec = std::from_chars(c, last, data[1]);
  if (pec.ec != std::errc{}) {
    return 1;
  }

  // skip one column, aka go to next comma
  c = pec.ptr;
  while (*++c && *c != ',')
    ;
  if (!*c) {
    return 1;
  }
  pec = std::from_chars(++c, last, data[2]);
  if (pec.ec != std::errc{})
    return 1;

  c = pec.ptr;
  if (*c++ != ',')
    return 1;
  pec = std::from_chars(c, last, data[3]);
  if (pec.ec != std::errc{})
    return 1;

  c = pec.ptr;
  if (*c++ != ',')
    return 1;
  pec = std::from_chars(c, last, data[4]);
  if (pec.ec != std::errc{})
    return 1;

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
    dso::utils::celestrak::details::CelestTrakSWFlux &flux_data) noexcept {

  // open the csv file
  std::ifstream fin(fncsv);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR. Failed opening SW csv file %s\n", fncsv);
    return fn_error;
  }

  // read first line; should match the expected header
  char line[MAX_SW_CHARS];
  const auto hsz = std::strlen(header);
  if (!fin.getline(line, MAX_SW_CHARS)) {
    return fn_error;
    // may include whitespace chars ... no std::strcmp
  } else if (std::strncmp(line, header, hsz)) {
    fprintf(stderr, "ERROR. Failed to match header line in SW csv file %s\n",
            fncsv);
    return fn_error;
  }

  dso::modified_julian_day target = mjd;
  dso::modified_julian_day cmjd;

  // read first data line; if record date > target data, the target date is
  // not included in the file ...
  fin.getline(line, MAX_SW_CHARS);
  if (resolve_csv_line_date(line, cmjd)) {
    fprintf(stderr, "ERROR. Failed parsing SW csv line %s\n", line);
    return parse_error;
  }
  if (cmjd > target)
    return target_prior;

  // keep on parsing/reading lines, untill we reach target date or EOF
  int error = 0;
  while (cmjd != target && !error) {
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

  // date matched! relevant record is stored in line
  double data[5] = {0};
  if (resolve_csv_line_records(line, data)) {
    fprintf(stderr, "ERROR. Failed resolving flux data from line %s\n", line);
    return parse_error;
  }

  flux_data.mjd_ = cmjd;
  flux_data.f107Obs = data[0];
  flux_data.f107Adj = data[1];
  flux_data.f107ObsC81 = data[2];
  flux_data.f107ObsL81 = data[3];
  flux_data.f107AdjC81 = data[4];

  return 0;
}
