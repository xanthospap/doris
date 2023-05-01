#include "tides.hpp"
#include <algorithm>
#include <charconv>
#include <cstring>
#include <exception>
#include <fstream>
#include <limits>
#include <vector>

constexpr const int MAX_CHARS_IN_OCCOEFFS = 256;

namespace {
struct OcTideRecordLine {
  /* Doodson number */
  dso::DoodsonNumber doodson;
  /* l: degree m: order */
  int l, m;
  /* coefficients: C_lm^(+), S_lm^(+), C_lm^(-), S_lm^(-)
   * i.e.          DelC+     DelS+     DelC-     DelS-
   */
  double DelCpl, DelSpl, DelCmi, DelSmi;
}; /* OcTideRecordLine */

/* @brief Parse a Doodson is/number from a line of type:
 * Doodson Darw  l   m    DelC+     DelS+       DelC-     DelS-
 * 055.565 om1   1   0  -0.84987   0.00000    -0.84987   0.00000
 * The function reads and examines only the first 7 characters. The rest of
 * the line is not considered.
 * @param[in] line A line where the first 7 characters hold a Doodson number
 * @param[out] d   The resolved Doodson number as dso::DoodsonNumber, if the
 *                 returned status is ok
 */
dso::iStatus parse_doodson_iers_str(const char *str,
                                    dso::DoodsonNumber &d) noexcept {
  /* check that the input string has size >=7 and that all but the 3rd char
   * are actually numeric values
   */
  for (int i = 0; i < 7; i++) {
    if (!str[i] || !(std::isdigit(*(unsigned char *)(str + i)) || i == 3)) {
      fprintf(stderr,
              "[ERROR] Failed to resolve Doodson number from string %s "
              "(traceback: %s)\n",
              str, __func__);
      return dso::iStatus(1);
    }
  }

  /* first three ints, for variables τ, s, h */
  d(0) = str[0] - '0';
  d(1) = str[1] - '0' - 5;
  d(2) = str[2] - '0' - 5;

  /* 4th character should be either a '.' or a ',' */
  if (str[3] != '.' && str[3] != ',') {
    fprintf(stderr,
            "[ERROR] Failed to resolve Doodson number from string %s "
            "(traceback: %s)\n",
            str, __func__);
    return dso::iStatus(1);
  }

  /* next three characters, for variables p, N', ps */
  d(3) = str[4] - '0' - 5;
  d(4) = str[5] - '0' - 5;
  d(5) = str[6] - '0' - 5;

  /* all done */
  return dso::iStatus(0);
}

/* @brief Parse a OcTideRecordLine instance off from a line of type:
 * Doodson Darw  l   m    DelC+     DelS+       DelC-     DelS-
 * 055.565 om1   1   0  -0.84987   0.00000    -0.84987   0.00000
 * @param[in] line A line as above
 * @param[out] rec The resolved OcTideRecordLine recorded in the line, if the
 *                 returned status is ok
 */
dso::iStatus parse_line(const char *line, OcTideRecordLine &rec) noexcept {
  /* lambda function: skip whitespace characters */
  auto skipws = [](const char *str) noexcept -> const char * {
    while (*str && *str == ' ')
      ++str;
    return str;
  };

  /* read and resolve Doodson number. IERS-2010 uses +5 for the Doodson
   * number elements. these are the first 7 chars
   */
  if (parse_doodson_iers_str(skipws(line), rec.doodson)) {
    fprintf(stderr,
            "[ERROR] Failed to resolve Doodson number from line: %s "
            "(traceback: %s)\n",
            line, __func__);
    return dso::iStatus(1);
  }

  /* size of line */
  const int sz = std::strlen(line);

  /* error status */
  int error = 0;

  /* start at: 7(=Doodson Number) + 1(whitespace) */
  const char *str = line + 7 + 1;
  /* next four chars are the Darwin name; just ignore for now! */
  str += 4;

  /* read degree and order (two ints) */
  {
    auto cres = std::from_chars(skipws(str), line + sz, rec.l);
    error += (cres.ec != std::errc{});
    str = cres.ptr;
    cres = std::from_chars(skipws(str), line + sz, rec.m);
    error += (cres.ec != std::errc{});
    str = cres.ptr;
    if (error || (rec.l < 0 || rec.m > rec.l)) {
      fprintf(stderr,
              "[ERROR] Failed to resolve SH degree/order from line: %s "
              "(traceback: %s)\n",
              line, __func__);
      return dso::iStatus(1);
    }
  }

  /* read coefficients in the order: DelC+ DelS+ DelC- DelS- */
  double d[4];
  {
    for (int i = 0; i < 4; i++) {
      auto cres = std::from_chars(skipws(str), line + sz, d[i]);
      error += (cres.ec != std::errc{});
      str = cres.ptr;
    }
    if (error) {
      fprintf(stderr,
              "[ERROR] Failed to resolve SH coefficients from line: %s "
              "(traceback: %s)\n",
              line, __func__);
      return dso::iStatus(1);
    }
  }

  /* assign to passed-in instance */
  rec.DelCpl = d[0];
  rec.DelSpl = d[1];
  rec.DelCmi = d[2];
  rec.DelSmi = d[3];

  return dso::iStatus(0);
}

/// @note The input argument vec cannot be declared const, cause then the
///       return type of the find_if call will be an iterator to const,
///       hence different from guess
auto find_doodson_in_vec(
    std::vector<dso::DoodsonOceanTideConstituent> &vec,
    const dso::DoodsonNumber &d,
    std::vector<dso::DoodsonOceanTideConstituent>::iterator guess) noexcept {

  if (guess != vec.end() && guess->doodson_number() == d) [[likely]]
    return guess;

  auto it = std::find_if(
      vec.begin(), vec.end(),
      [D = std::cref(d)](const dso::DoodsonOceanTideConstituent &element) {
        return (element.doodson_number() == D);
      });
  /*[=](const dso::DoodsonOceanTideConstituent &el) {
    return el.doodson_number() == d;
  });*/

  return it;
}

} /* unnamed namespace */

/* @brief Read in and parse a IERS-distributed Cnm/Snm (+/-) coefficients for
 *        and ocean tide model. Example file:
 *        IERS2010/chapter6/tidemodels/fes2004_Cnm-Snm.dat
 *  This file contains coefficients:
 *  Doodson Darw  l   m    DelC+     DelS+       DelC-     DelS-
 *   55.565 Om1   2   0  -6.58128   0.00000    -0.00000  -0.00000
 *   [ ... ]
 *
 * @param[in] fn     Filename of the file to be parsed
 * @param[out] freqs A vector of DoodsonOceanTideConstituent instances, one
 *                   for each wave/constituent included in the input file
 * @param[in] max_degree Max Cnm/Snm degree to be collected (inclusive)
 * @param[in] max_order  Max Cnm/Snm order to be collected (inclusive)
 * @param[in] scale   Scale factor for the Cnm/Snm coefficients
 */
dso::iStatus dso::memmap_octide_coefficients(
    const char *fn, std::vector<dso::DoodsonOceanTideConstituent> &freqs,
    int max_degree, int max_order, double scale) noexcept {
  /* open input file */
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Failed opening file: %s (traceback: %s)\n", fn,
            __func__);
    return dso::iStatus(1);
  }

  freqs.clear();
  char line[MAX_CHARS_IN_OCCOEFFS];
  const char *HEADER =
      "Doodson Darw  l   m    DelC+     DelS+       DelC-     DelS-";

  /* Usually first 3 lines are comments and the next is header; skip as many
   * lines needed, to match the header
   */
  {
    int dummy_it = 0;
    int max_dummy_it = 15;
    while (++dummy_it < max_dummy_it && std::strcmp(line, HEADER))
      fin.getline(line, MAX_CHARS_IN_OCCOEFFS);
    if (dummy_it >= max_dummy_it) {
      fprintf(
          stderr,
          "[ERROR] Failed matching header line in file %s (traceback: %s)\n",
          fn, __func__);
      return dso::iStatus(1);
    }
  }

  /* i want to hold max degree and order that we can check it later on.
   * we will assume that all waves have the same max degree and order
   */
  int maxDegreeRead = std::numeric_limits<double>::min();
  int maxOrderRead = std::numeric_limits<double>::min();

  /* keep on reading constituents ... */
  OcTideRecordLine rec;
  auto it = freqs.end();
  dso::iStatus error = dso::iStatus(0);
  while (fin.getline(line, MAX_CHARS_IN_OCCOEFFS) && !error) {
    /* parse line holdings to rec */
    error = parse_line(line, rec);
    if (!error) {
      /* find constintuent in vector */
      it = find_doodson_in_vec(freqs, rec.doodson, it);
      if (it != freqs.cend()) {
        /* constituent found, add elements for given degree and order
         * (if needed)
         */
        if (rec.l <= max_degree && rec.m <= max_order) {
          it->delCp(rec.l, rec.m) = rec.DelCpl * scale;
          it->delSp(rec.l, rec.m) = rec.DelSpl * scale;
          it->delCm(rec.l, rec.m) = rec.DelCmi * scale;
          it->delSm(rec.l, rec.m) = rec.DelSmi * scale;
        }
      } else {
        /* else, add new constituent */
        freqs.emplace_back(rec.doodson, max_degree, max_order);
        it = std::prev(freqs.end());
        /* clear (i.e. set to zero) the new wave's Cnm/Snm coefs */
        it->clear_coefficients();
        if (rec.l <= max_degree && rec.m <= max_order) {
          it->delCp(rec.l, rec.m) = rec.DelCpl * scale;
          it->delSp(rec.l, rec.m) = rec.DelSpl * scale;
          it->delCm(rec.l, rec.m) = rec.DelCmi * scale;
          it->delSm(rec.l, rec.m) = rec.DelSmi * scale;
        }
      }
      if (rec.l >= maxDegreeRead)
        maxDegreeRead = rec.l;
      if (rec.m >= maxOrderRead)
        maxOrderRead = rec.m;
    } else {
      /* some kind of error while paring the line! */
      fprintf(stderr,
              "[ERROR] Failed parsing line in file %s (traceback: %s)\n", fn,
              __func__);
      fprintf(stderr, "[     ] line: %s\n", line);
      return dso::iStatus(1);
    }
  } /* end reading lines */

  /* let's check the degree and order of the waves collected */
  if (maxDegreeRead != max_degree) {
    fprintf(stderr,
            "[ERROR] Requested paring constituents with max degree = %d but "
            "only found coefficients up to degree %d (traceback: %s)\n",
            max_degree, maxDegreeRead, __func__);
    return dso::iStatus(1);
  }
  if (maxOrderRead != max_order) {
    fprintf(stderr,
            "[ERROR] Requested paring constituents with max order = %d but "
            "only found coefficients up to order %d (traceback: %s)\n",
            max_order, maxOrderRead, __func__);
    return dso::iStatus(1);
  }

  /* have we read the file through ? */
  if (!fin.good() && fin.eof())
    return dso::iStatus(0);

  return dso::iStatus(1);
}
