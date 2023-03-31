#include "atmosphere/dtm2020/dtm2020.hpp"
#include <charconv>
#include <cstdio>
#include <fstream>
#include <cstring>

namespace {
const char *skip_ws(const char *str) noexcept {
  while (*str && std::isspace(*str))
    ++str;
  return str;
}

/* @brief Parse numeric values off from a file fescribed in
 * https://github.com/swami-h2020-eu/mcm/blob/main/data/DTM_2020_F107_Kp.dat
 * Resolved values in order:
 * TT H HE O N2 O2 N  T0 TP are stored in the data array in indexes [0-9]
 *
 * @param[in] line A (data) line from the file 'DTM_2020_F107_Kp.dat' but
 *                 not the first two (which are header)
 * @param[in] data An array of size >= 9 to store line values
 * @return A non-zero return value, denotes an error.
 */
int parse_data_line(const char *line, double *data) noexcept {
  const int sz = std::strlen(line);
  int num;
  double sigma;

  /* Format of data file is:
   * https://github.com/swami-h2020-eu/mcm/blob/main/data/DTM_2020_F107_Kp.dat
   * num TT H HE O N2 O2 N  T0 TP
   * 0   1  3 5  7 9  11 13 15 17
   * each value is followed by its sigma value, which is not of interest
   * columns are: 1(a/a) + 9 * 2 (values+sigmas)
   */
  auto cres = std::from_chars(skip_ws(line), line+sz, num); /* discarded */
  if (cres.ec != std::errc{}) {
    fprintf(stderr,
            "[ERROR] Failed reading line number in DTM2020 data file "
            "(traceback: %s)\n",
            __func__);
    return 1;
  }

  int error = 0;
  for (int i = 0; i < 9 && (!error); i++) {
    line = cres.ptr;
    cres = std::from_chars(skip_ws(line), line+sz, data[i]); /* read value */
    error += (cres.ec != std::errc{});
    line = cres.ptr;
    cres = std::from_chars(skip_ws(line), line+sz, sigma); /* read sigma */
    error += (cres.ec != std::errc{});
  }

  if (error) {
    fprintf(stderr,
            "[ERROR] Failed resolving data line from DTM2020 data file "
            "(traceback: %s)\n",
            __func__);
    return 1;
  }
  return 0;
}
} // unnamed namespace

/*
 * This is a translated version of the FORTAN source code found at
 * https://github.com/swami-h2020-eu/mcm/blob/main/src/dtm2020/dtm2020_F107_Kp-subr.f90
 * Original (FORTRAN) Author and translated Version:
 * Author : SB (CNES)
 * last version updated : 21/10/2020
 *
 * This function reads in and parses the fields:
 * TT H HE O N2 O2 N  T0 TP
 * in the respective arrays of the calling instance (az,az2,h,he,o,o2,t0,tp,tt)
 * 
 * @param[in] fn The data file DTM_2020_F107_Kp.dat
 */
int dso::Dtm2020::map_model_coeffs(const char *fn) noexcept {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(
        stderr,
        "[ERROR] Failed opening DTM2020 coefficient file %s (traceback: %s)\n",
        fn, __func__);
    return 1;
  }

  constexpr const int MAX_CHARS = 512;
  char line[MAX_CHARS];

  /* status */
  int error = 0;

  /* first line is header, just skip it (don't even check) */
  fin.getline(line, MAX_CHARS);
  const char *end = line + std::strlen(line);

  /* second line includes (first numeric) the number of records */
  int num_lines = 0;
  fin.getline(line, MAX_CHARS);
  auto cres = std::from_chars(skip_ws(line), end, num_lines);
  if (cres.ec != std::errc{} || error || num_lines <= 0) {
    fprintf(stderr,
            "[ERROR] Failed resolving number of data lines in DTM2020 "
            "coefficient file %s (traceback: %s)\n",
            fn, __func__);
    return 1;
  }

  /* well, num_lines should be at maximum nlatm */
  assert(num_lines <= nlatm);

  /* so ... keep on reading and resolving ... */
  double data[9];
  for (int i = 0; i < num_lines && (!error); i++) {
    fin.getline(line, MAX_CHARS);
    error = parse_data_line(line, data);
    /* tt(i),h(i),he(i),o(i),az2(i),o2(i),az(i),t0(i),tp(i) */
    tt[i] = data[0];
    h[i] = data[1];
    he[i] = data[2];
    o[i] = data[3];
    az2[i] = data[4];
    o2[i] = data[5];
    az[i] = data[6];
    t0[i] = data[7];
    tp[i] = data[8];
  }

  if (error) {
    fprintf(stderr, "[ERROR] Failed parsing data from file %s (taceback: %s)\n",
            fn, __func__);
    return 1;
  }
  return 0;
}
