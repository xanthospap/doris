#include "datetime/dtcalendar.hpp"
#include <datetime/dtfund.hpp>
#include <fstream>
#include <cstdio>
#include <system_error>
#include <charconv>
#include <limits>
#include <cstring>

/*
https://celestrak.org/SpaceData/
https://pypi.org/project/nrlmsise00/
*/

constexpr const double MissingSwData = std::numeric_limits<double>::min();
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

struct CsvLineInfo {
  dso::modified_julian_day mjd_;
  double f107Obs{MissingSwData}, f107Adj{MissingSwData},
      f107ObsC81{MissingSwData}, f107ObsL81{MissingSwData},
      f107AdjC81{MissingSwData};
};

int resolve_csv_line_records(const char *line, double *data) noexcept {
  const auto sz = std::strlen(line);

  printf("calling %s\n", __func__);
  printf("resolving line %s", line);
  return 444;
  
  // find and skip 24 ',' chars
  int commas = 0;
  const char *c = line;
  while (*c && commas<24) {
    if (*c++ == ',') {
      ++commas;
    }
  }
  if (!*c) return 1;
  printf("cool!\n");
  return 555;

  // find next comma offset
  //int offset = 0;
  //const char *end = c;
  //while (end && *end++ != ',') {}
  //if (!end) return 1;

  // 
  const char *last = line + sz - 1;

  printf("parsing line: [%s", c);

  // resolve ...
  auto pec = std::from_chars(c, last, data[0]);
  if (pec.ec != std::errc{})
    return 1;

  // now ptr should point to the first, non-resolved char, aka next comma!
  c = pec.ptr;
  if (*c != ',') return 1;
  pec = std::from_chars(c, last, data[1]);
  if (pec.ec != std::errc{})
    return 1;

  while (*c && *c++ != ',');
  if (!*c) return 1;
  pec = std::from_chars(c, last, data[2]);
  if (pec.ec != std::errc{})
    return 1;
  
  c = pec.ptr;
  if (*c != ',') return 1;
  pec = std::from_chars(c, last, data[3]);
  if (pec.ec != std::errc{})
    return 1;

  c = pec.ptr;
  if (*c != ',') return 1;
  pec = std::from_chars(c, last, data[4]);
  if (pec.ec != std::errc{})
    return 1;

  return 0;
}


int resolve_csv_line_date(const char *line,
                          dso::modified_julian_day &mjd) noexcept {
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

int parse_csv_for_date(const dso::datetime<dso::seconds> &t,
                       const char *fncsv, CsvLineInfo &flux_data) noexcept {

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

  dso::modified_julian_day target = t.mjd();
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
    printf("parsed line %s\n", line);
  }

  // check for error/EOF
  if (error) {
    if (fin.eof()) return eof_reached;
    return error;
  }

  // date matched! relevant record is stored in line
  double data[5];
  if (resolve_csv_line_records(line, data)) {
    fprintf(stderr, "ERROR. Failed resolving flux data from line %s", line);
    return parse_error;
  }

  flux_data.mjd_ = cmjd;
  flux_data.f107Obs   = data[0];
  flux_data.f107Adj   = data[1];
  flux_data.f107ObsC81= data[2];
  flux_data.f107ObsL81= data[3];
  flux_data.f107AdjC81= data[4];

  return 0;
}

int main(int argc, char *argv[]) {
  if (argc!=2) {
    fprintf(stderr, "Usage %s [CSV FLUX FILE]\n", argv[0]);
    return 1;
  }

  dso::datetime<dso::seconds> tt(dso::year(2022), dso::month(1),
                                 dso::day_of_month(1), dso::seconds());

  CsvLineInfo flux_data;
  if (int error; (error=parse_csv_for_date(tt, argv[1], flux_data))) {
    //fprintf(stderr, "ERROR. Failed with error code: %d\n", error);
    return error;
  }

  printf("Flux data for requested date:\n");
  printf("F10.7_OBS         = %+.2f\n", flux_data.f107Obs);
  printf("F10.7_ADJ         = %+.2f\n", flux_data.f107Adj);
  printf("F10.7_OBS_CENTER81= %+.2f\n", flux_data.f107ObsC81);
  printf("F10.7_OBS_LAST81  = %+.2f\n", flux_data.f107ObsL81);
  printf("F10.7_ADJ_CENTER81= %+.2f\n", flux_data.f107AdjC81);

  return 0;
}
