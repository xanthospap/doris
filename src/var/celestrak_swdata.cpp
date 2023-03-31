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

int dso::utils::celestrak::details::resolve_csv_line_records(
    const char *line,
    dso::utils::celestrak::details::CelestTrakSWFlux &flux_data) noexcept {

  /* parse date UTC */
  if (dso::utils::celestrak::details::resolve_csv_line_date(line,
                                                            flux_data.mjd_)) {
    return 1;
  }

  const auto sz = std::strlen(line);
  const char *c = line;

  /* find and skip 3 ',' chars i.e. columns [DATE,BSRN,ND,] */
  int commas = 0;
  while (*c++ && commas < 3) {
    if (*c == ',') {
      ++commas;
    }
  }
  if (!*c) {
    return 2;
  }

  /* last character in the string, the limit for conversions */
  const char *last = line + sz - 1;

  /* get 8 hourly Kp's (nts) i.e. [KP1,KP2,KP3,KP4,KP5,KP6,KP7,KP8,] */
  int i=0;
  while (i<8) {
    auto pec = std::from_chars(c, last, flux_data.KpIndexes[i]);
    c = pec.ptr;
    if (*c++ != ',' || pec.ec != std::errc{}) {
      return 3 + i;
    }
    ++i;
  }

  /* get 8 hourly Kp's (nts) i.e. [KPSUM,] */
    auto pec = std::from_chars(c, last, flux_data.KpSum);
    c = pec.ptr;
    if (*c++ != ',' || pec.ec != std::errc{}) {
      return 9;
    }

  /* get 8 hourly Ap's (nts) i.e. [AP1,AP2,AP3,AP4,AP5,AP6,AP7,AP8,] */
  while (i < 8) {
    pec = std::from_chars(c, last, flux_data.ApIndexes[i]);
    c = pec.ptr;
    if (*c++ != ',' || pec.ec != std::errc{}) {
      return 3 + i;
    }
    ++i;
  }

  /* get Ap daily average [AP_AVG,] */
    pec = std::from_chars(c, last, flux_data.ApDailyAverage);
    c = pec.ptr;
    if (*c++ != ',' || pec.ec != std::errc{}) {
      return 11;
    }

  {/* skip columns [CP,C9,ISN,] */
    commas = 0;
    while (*c++ && commas < 3) {
      if (*c == ',') {
        ++commas;
      }
    }
    if (!*c)
      return 12;
  }

  /* we are not at the F10.7_OBS column */
    pec = std::from_chars(c, last, flux_data.f107Obs);
    c = pec.ptr;
    if (pec.ec != std::errc{} || *c++ != ',') {
      return 13;
    }
  
  /* get [F10.7_ADJ,] */
  pec = std::from_chars(c, last, flux_data.f107Adj);
  c = pec.ptr;
  if (pec.ec != std::errc{} || *c++ != ',') {
    return 14;
  }

  /* [F10.7_DATA_TYPE] just copy three chars */
  std::memcpy(flux_data.flag, c, 3 * sizeof(char));
  c += 3;
  if (*c++ != ',')
    return 15;

  /* [F10.7_OBS_CENTER81,] */
  pec = std::from_chars(c, last, flux_data.f107ObsC81);
  c = pec.ptr;
  if (pec.ec != std::errc{} || *c++ != ',')
    return 16;

  /* [F10.7_OBS_LAST81,] */
  pec = std::from_chars(c, last, flux_data.f107ObsL81);
  c = pec.ptr;
  if (pec.ec != std::errc{} || *c++ != ',')
    return 17;

  /* [F10.7_ADJ_CENTER81,] */
  pec = std::from_chars(c, last, flux_data.f107AdjC81);
  c = pec.ptr;
  if (pec.ec != std::errc{} || *c++ != ',')
    return 18;
  
  /* [F10.7_ADJ_LAST81] note that this is the last column!*/
  pec = std::from_chars(c, last, flux_data.f107AdjL81);
  c = pec.ptr;
  if (pec.ec != std::errc{})
    return 18;

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

  // to datetime ... (note that this is UTC)
  const dso::datetime<dso::seconds> d(dso::year(year), dso::month(month),
                                      dso::day_of_month(day), dso::seconds(0));
  mjd = d.mjd();

  // all done
  return 0;
}

int dso::utils::celestrak::details::get_next(
    dso::modified_julian_day mjd, const char *fncsv,
    CelestTrakSWFlux *flux_data, int flux_data_sz,
    std::ifstream::pos_type &fpos) noexcept {
  // try getting the data for the given mjd; store it in new_flux
  // fpos set to new, latest position in file
  dso::utils::celestrak::details::CelestTrakSWFlux new_flux;
  if (parse_csv_for_date(mjd, fncsv, &new_flux, fpos, 0, 0)) {
    return 1;
  }

  CelestTrakSWFlux *__restrict__ fp = flux_data;
  auto fsz = sizeof(CelestTrakSWFlux);

  // left-shit elements in array
  for (int i = 1; i < flux_data_sz; i++) {
    std::memcpy(fp - 1 + i, fp + i, fsz);
  }

  // store the retrieved data in the rightmost element
  std::memcpy(fp + flux_data_sz - 1, &new_flux, fsz);

#ifdef DEBUG
  for (int i = 0; i < flux_data_sz - 1; i++) {
    assert(flux_data[i].mjd_ + dso::modified_julian_day(1) ==
           flux_data[i + 1].mjd_);
  }
#endif

  // all done
  return 0;
}

/* TODO obsolete, should be removed */
int dso::utils::celestrak::details::parse_csv_for_date(
    dso::modified_julian_day mjd, const char *fncsv,
    dso::utils::celestrak::details::CelestTrakSWFlux *flux_data,
    std::ifstream::pos_type &fpos, int days_before, int days_after) noexcept {

  // set the fpos at the start of the file and save the variable in fpos2
  const auto fpos2 = fpos;
  // !only set this to anything other than 0, if the dates are resolved
  // successefully!
  fpos = 0;

  // open the csv file
  std::ifstream fin(fncsv);
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR]: Failed opening SW csv file %s (traceback: %s)\n",
            fncsv, __func__);
    return fn_error;
  }
  char line[MAX_SW_CHARS];

  // set mjd's
  dso::modified_julian_day start_mjd(mjd.as_underlying_type() - days_before);
  dso::modified_julian_day end_mjd(mjd.as_underlying_type() + days_after);

  // ok, re-set fpos to where it was, ang go there
  fin.seekg(fpos2);
  if (!fpos2) {
    // we are reading from the top! read first line; should match the expected
    // header
    const auto hsz = std::strlen(header);
    if (!fin.getline(line, MAX_SW_CHARS)) {
      return fn_error;
      // may include whitespace chars ... no std::strcmp
    } else if (std::strncmp(line, header, hsz)) {
      fprintf(stderr,
              "[ERROR] Failed to match header line in SW csv file %s (traceback: %s)\n",
              fncsv, __func__);
      return fn_error;
    }
  }

  dso::modified_julian_day cmjd;

  // read first data line; if record_date > start_mjd data, the target date is
  // not included in the file ...
  fin.getline(line, MAX_SW_CHARS);
  if (resolve_csv_line_date(line, cmjd)) {
    fprintf(stderr, "[ERROR]: Failed parsing SW csv line %s (traceback: %s)\n", __func__,
            line);
    return parse_error;
  }
  if (cmjd > start_mjd) {
    // wait! don't give up! if fpos was somethin other than 0, maybe it was
    // erronuous, and the date exists in the file but is before fpos. try the
    // function once more, with an fpos=0
    if (fpos2)
      return parse_csv_for_date(mjd, fncsv, flux_data, fpos, days_before,
                                days_after);
    return target_prior;
  }

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

  // set fpos
  fpos = fin.tellg();

  return 0;
}
