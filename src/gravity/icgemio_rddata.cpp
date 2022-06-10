#include "icgemio.hpp"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

///< approximate max data line length, see
///< http://icgem.gfz-potsdam.de/ICGEM-Format-2011.pdf
constexpr std::size_t max_data_line = 512;

bool starts_with(const char *pattern, const char *line) noexcept {
  if (line)
    return !std::strncmp(pattern, line, std::strlen(pattern));
  return false;
}

std::size_t coeffs_nr(int l, int m) noexcept {
  if (l == m) {
    std::size_t n = l + 1;
    return (n * (n - 1)) / 2 + n;
  }

  std::size_t sum = 0;
  for (int i = 0; i <= l; i++)
    sum += (i > m) ? (m + 1) : (i + 1);
  return sum;
}

/// @warning coeffs should have already been initialized and allocated with
///          enough memmory to hold the (to-be-) parsed coefficients.
int dso::Icgem::parse_data(int l, int m, dso::HarmonicCoeffs *coeffs) noexcept {

  int error = 0;
  if (l > max_degree || m > l) {
    fprintf(
        stderr,
        "[ERROR] Invalid degree/order given to data parse (traceback: %s)\n",
        __func__);
    error = 1;
  }

  if (coeffs->degree() < l) {
    fprintf(stderr,
            "[ERROR] Cannot read harmonics of degree %d to HarmonicsCoeffs of "
            "degree %d (traceback: %s)\n",
            l, coeffs->degree(), __func__);
    error = 1;
  }

  std::ifstream fin(filename.c_str());
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Failed opening icgem file %s (traceback: %s)\n",
            filename.c_str(), __func__);
    error = 1;
  }

  if (error) return error;

  // assign gravity model constants
  coeffs->GM() = earth_gravity_constant;
  coeffs->Re() = radius;
  coeffs->normalized() = this->is_normalized();

  fin.seekg(data_section_pos);

  // Note (1)
  // -------------------------------------------------------------------------
  // for some gfc files (e.g. the EGM2008) it may happen that the values for
  // C_10 and C_11 are missing; that is because the are nominally zero; here
  // we will set them to some random value, so that if at the end of parsing
  // we are missing exactly two values and the C_10 and C_11 have the values
  // set here, these are the ones not parsed!
  coeffs->C(1, 0) = -999e0;
  coeffs->C(1, 1) = -999e0;

  char line[max_data_line], *end, *start;
  int ll, mm;
  double Clm, Slm;
  std::size_t coeffs_read = 0, coeffs_to_read = coeffs_nr(l, m);

  fin.getline(line, max_data_line);
  while (fin.good() && coeffs_read < coeffs_to_read) {
    // we are only interested in lines that start with 'gfc' (but not 'gfct')
    if (starts_with("gfc ", line)) {

      // expecting columns: degree, order, Clm, Slm, [...]; note that it
      // (seldom) happens that the doubles are written in fortran format ...
      start = line + 4;
      ll = std::strtol(start, &end, 10);
      if (end == start) {
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      start = end;
      mm = std::strtol(start, &end, 10);
      if (end == start) {
        fprintf(stderr,
                "[ERROR] Failed parsing order parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      // only interested in the coefficients, if degree and order are less than
      // max
      if (ll <= l && mm <= m) {

        start = end;
        Clm = std::strtod(start, &end);
        if (end == start) {
          fprintf(stderr,
                  "[ERROR] Failed parsing Clm parameter in line: [%s]; icgem "
                  "file %s (traceback: %s)\n",
                  line, filename.c_str(), __func__);
          return 1;
        }

        start = end;
        Slm = std::strtod(start, &end);
        if (end == start) {
          fprintf(stderr,
                  "[ERROR] Failed parsing Slm parameter in line: [%s]; icgem "
                  "file %s (traceback: %s)\n",
                  line, filename.c_str(), __func__);
          return 1;
        }

        // assign to harmonic coefficients matrix
        coeffs->C(ll, mm) = Clm;
        ++coeffs_read;
        if (mm == 0) {
          assert(Slm == 0e0);
        } else {
          coeffs->S(ll, mm) = Slm;
        }
      }
    }

    fin.getline(line, max_data_line);
  }

  if (coeffs_read < coeffs_to_read) {
    // before reporting an error, see if we are in the case described in
    // note (1)
    if (coeffs_to_read - coeffs_read == 2 &&
        (coeffs->C(1, 0) == coeffs->C(1, 1) && coeffs->C(1, 1) == -999e0)) {
      printf("[NOTE] The coefficients C(1,0) and C(1,1) are not explicitly "
             "written in the icgem file: %s\n",
             filename.c_str());
      printf("[NOTE] Setting C(1, 0) = C(1, 1) = 0e0\n");
      coeffs->C(1, 0) = 0e0;
      coeffs->C(1, 1) = 0e0;
    } else {
      fprintf(stderr,
              "[ERROR] EOF reached before reading all Snm/Cnm coefficients! "
              "read/expected %lu/%lu; icgem file %s (traceback: %s)\n",
              coeffs_read, coeffs_to_read, filename.c_str(), __func__);
      return 2;
    }
  }

  return 0;
}
