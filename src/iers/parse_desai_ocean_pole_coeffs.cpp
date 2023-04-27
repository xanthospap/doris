#include "cmat2d.hpp"
#include "base_error.hpp"
#include "tides.hpp"
#include <cctype>
#include <charconv>
#include <cstdio>
#include <fstream>

/* Parse ocean pole tide coefficients provided by IERS, see
 * IERS 2010, Sec. 6.4
 * ftp://tai.bipm.org/iers/conv2010/chapter6/desaiscopolecoef.txt
 */

namespace {
/* max chars in input file per line */
constexpr const int MAXSZ = 512;
/* header of the file (expected) missing whitespace characters */
const char *header = "nmAnm(Real)Bnm(Real)Anm(Imaginary)Bnm(Imaginary)";
/* skip space characters */
const char *skip_space(const char *line) noexcept {
  while (*line && std::isspace(line[0]))
    ++line;
  return line;
}
} // namespace

dso::iStatus dso::parse_desai_ocean_pole_tide_coeffs(
    const char *fn, int max_degree, int max_order,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &Areal,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &Aimag,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &Breal,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &Bimag) noexcept {
  /* fill output arrays with zeros */
  Areal.fill_with(0e0);
  Aimag.fill_with(0e0);
  Breal.fill_with(0e0);
  Bimag.fill_with(0e0);

  std::fstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Failed to open ocean pole tide coefficients file %s "
            "(traceback: %s)\n",
            fn, __func__);
    return dso::iStatus(1);
  }

  /* size validation */
  if (max_degree > 360 || max_order > 360 || max_order > max_degree ||
      Areal.rows() < max_degree || Areal.cols() < max_order) {
    fprintf(stderr,
            "[ERROR] Erronuous input sizes given for parsing ocean pole tide "
            "coefficients "
            "(traceback: %s)\n",
            __func__);
    return dso::iStatus(1);
  }

  /* line buffer */
  char line[MAXSZ];

  /* check header */
  {
    fin.getline(line, MAXSZ);
    /* copy the line to hdr ommiting spaces */
    char hdr[MAXSZ];
    const char *c = line;
    int i = 0;
    while (*c) {
      if (!std::isspace(*c)) {
        hdr[i] = *c;
        ++i;
      }
      ++c;
    }
    hdr[i] = '\0';
    /* check if space-less header is ok */
    if (std::strcmp(hdr, header)) {
      fprintf(stderr,
              "[ERROR] Failed to validate header in ocean pole tide "
              "coefficients file %s (traceback: %s)\n",
              fn, __func__);
      fprintf(stderr, "[ERROR] Read header line: [%s] (traceback: %s)\n", line,
              __func__);
      return dso::iStatus(1);
    }
  }

  /* go on, parsing data one line at a time */
  int n, m;
  double data[4];
  int error = 0;
  while (fin.getline(line, MAXSZ)) {
    const char *end = line + std::strlen(line);
    const char *c = line;
    /* parse current degree and order */
    auto res = std::from_chars(skip_space(c), end, n);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    res = std::from_chars(skip_space(c), end, m);
    if (res.ec != std::errc{})
      ++error;
    c = res.ptr;
    if (error) {
      fprintf(stderr,
              "[ERROR] Failed parsing data from ocean pole tide coefficients "
              "file %s (traceback: %s)\n",
              fn, __func__);
      fprintf(stderr, "[ERROR] line was [%s] (traceback: %s)\n", line,
              __func__);
      return dso::iStatus(1);
    }
    /* if needed parse coefficients */
    if (n <= max_degree && m <= max_order) {
      for (int i = 0; i < 4; i++) {
        res = std::from_chars(skip_space(c), end, data[i]);
        if (res.ec != std::errc{})
          ++error;
        c = res.ptr;
      }
      if (error) {
        fprintf(stderr,
                "[ERROR] Failed parsing data from ocean pole tide coefficients "
                "file %s (traceback: %s)\n",
                fn, __func__);
        fprintf(stderr, "[ERROR] line was [%s] (traceback: %s)\n", line,
                __func__);
        return dso::iStatus(1);
      }
      /* assign */
      Areal(n, m) = data[0];
      Breal(n, m) = data[1];
      Aimag(n, m) = data[2];
      Bimag(n, m) = data[3];
    }
  } /* end looping through lines */

  /* all done */
  return dso::iStatus(0);
}
