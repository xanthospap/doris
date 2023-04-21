#include "solar_activity_data.hpp"
#include <cstdio>
#include <cstring>

using dso::utils::celestrak::details::resolve_csv_line_date;
using dso::utils::celestrak::details::resolve_csv_line_records;

namespace {
constexpr const int MAX_SW_CHARS = 512;
constexpr const char header[] =
    "DATE,BSRN,ND,KP1,KP2,KP3,KP4,KP5,KP6,KP7,KP8,KP_SUM,AP1,AP2,AP3,AP4,AP5,"
    "AP6,AP7,AP8,AP_AVG,CP,C9,ISN,F10.7_OBS,F10.7_ADJ,F10.7_DATA_TYPE,F10.7_"
    "OBS_CENTER81,F10.7_OBS_LAST81,F10.7_ADJ_CENTER81,F10.7_ADJ_LAST81";
} // unnamed namespace

int dso::SolarActivityData::feed(const char *fncsv, const dso::TwoPartDate &t,
                                 int days_before, int days_after) noexcept {
  /* clear output vector and reserve memmory */
  vdata.clear();
  vdata.reserve(days_before + days_after + 1);

  /* get central date's mjd */
  dso::modified_julian_day mjd(t.big());

  /* open the csv file */
  std::ifstream fin(fncsv);
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR]: Failed opening SW csv file %s (traceback: %s)\n",
            fncsv, __func__);
    return 1;
  }
  char line[MAX_SW_CHARS];

  /* set mjd's, collect data for [start_mjd, end_mjd) */
  dso::modified_julian_day start_mjd(mjd.as_underlying_type() - days_before);
  dso::modified_julian_day end_mjd(mjd.as_underlying_type() + days_after);

  /* should match the expected header */
  const auto hsz = std::strlen(header);
  if (!fin.getline(line, MAX_SW_CHARS)) {
    return 1;
    /* may include whitespace chars ... no std::strcmp */
  } else if (std::strncmp(line, header, hsz)) {
    fprintf(stderr,
            "[ERROR] Failed to match header line in SW csv file %s (traceback: "
            "%s)\n",
            fncsv, __func__);
    return 1;
  }

  /* current (i.e. while reading) MJD */
  dso::modified_julian_day cmjd;
  /* current flux data */
  dso::utils::celestrak::details::CelestTrakSWFlux fdata;
  /* keep on reading data unitll EOF or requested interval buffered */
  while (fin.getline(line, MAX_SW_CHARS)) {
    /* resolve date stored in line */
    if (int error = resolve_csv_line_date(line, cmjd); error) {
      fprintf(stderr,
              "[ERROR] Failed to resolve Solar Fulx data record from line %s "
              "(error=%d, traceback: %s)\n",
              line, error, __func__);
      return 1;
    }

    if (cmjd < start_mjd) {
      ;
    } else if (cmjd >= end_mjd) {
      break;
    } else {
      /* resolve data, store in flux_data[cdi] */
      if (int error = resolve_csv_line_records(line, fdata); error) {
        fprintf(stderr,
                "[ERROR] Failed resolving flux data from line %s (error=%d, "
                "traceback: %s)\n",
                line, error, __func__);
        return 1;
      }
      vdata.emplace_back(fdata);
    }
  }

  /* check for EOF */
  if (fin.eof()) {
    fprintf(stderr,
            "[ERROR] EOF encountered before finished reading solar flux data "
            "from %s (traceback: %s)\n",
            fncsv, __func__);
    return 2;
  }

  return 0;
}
