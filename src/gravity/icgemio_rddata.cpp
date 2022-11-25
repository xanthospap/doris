#include "icgemio.hpp"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <charconv>
#include <ggdatetime/dtfund.hpp>

///< approximate max data line length, see
///< http://icgem.gfz-potsdam.de/ICGEM-Format-2011.pdf
namespace {
constexpr std::size_t max_data_line = 512;

bool starts_with(const char *pattern, const char *line) noexcept {
  if (line)
    return !std::strncmp(pattern, line, std::strlen(pattern));
  return false;
}

const char *skip_ws(const char *str) noexcept {
  while (*str && *str==' ') ++str;
  return str;
}

int count_columns(const char *line) noexcept {
  const char *c = line;
  int cols = 0;
  while (*c) {
    while (*c && *c == ' ') ++c;
    if (*c && *c != ' ') ++cols;
    while (*c && *c != ' ') ++c;
  }
  return cols;
}

const char *goto_column(const char *line, int colnr) noexcept {
  const char *c = line;
  int cols = 0;
  while (*c) {
    while (*c && *c == ' ') ++c;
    if (*c && *c != ' ')
      if (++cols == colnr)
        return c;
    while (*c && *c != ' ') ++c;
  }
  return nullptr;
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

int resolve_date(const char *date_str,
                 dso::datetime<dso::nanoseconds> &t) noexcept {
  // expected format: t0[yyyymmdd.xxxx]
  using SecIntType = dso::nanoseconds::underlying_type;
  int year, month, dom, error = 0;
  double fraction, dummy;
  auto rc = std::from_chars(date_str, date_str + 4, year);
  error += (rc != std::errc{});
  rc = std::from_chars(date_str + 4, date_str + 6, month);
  error += (rc != std::errc{});
  rc = std::from_chars(date_str + 6, date_str + 8, dom);
  error += (rc != std::errc{});
  rc = std::from_chars(date_str, date_str + 14, fraction);
  fraction = std::modf(fraction, &dummy);
  const double fnanosec = dso::nanoseconds::sec_factor<double>();
  if (error) return error;
  t = dso::datetime<dso::nanoseconds>(
      dso::year(year), dso::month(month), dso::day_of_month(dom),
      dso::nanoseconds(static_cast<SecIntType>(fraction)));
  return 0;
}
} // unnamed namespace

/// @warning coeffs should have already been initialized and allocated with
///          enough memmory to hold the (to-be-) parsed coefficients.
int dso::Icgem::parse_data(int l, int m, const dso::datetime<dso::nanoseconds> &t, dso::HarmonicCoeffs *coeffs) noexcept {

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
  int error = 0;

  fin.getline(line, max_data_line);
  while (fin.good() && !error && coeffs_read < coeffs_to_read) {
    
    // gfc lines are for static-gravity field (if L=0=M, no effect ...)
    if (starts_with("gfc ", line)) {
      // expecting columns: degree, order, Clm, Slm, [...]; note that it
      // (seldom) happens that the doubles are written in fortran format ...
      start = line + 4;
      int ll;
      auto ccres = std::from_chars(skip_ws(start), line+sz, ll);
      if (ccres.ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        ++error;
      }

      int mm;
      ccres = std::from_chars(skip_ws(ccres.ptr), line+sz, mm);
      if (ccres.ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing order parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        ++error;
      }
      
      // only interested in the coefficients, if degree and order are less than
      // max
      if (ll <= l && mm <= m) {

        ccres = std::from_chars(skip_ws(ccres.ptr), line + sz, Clm);
        if (ccres.ec != std::errc{}) {
          fprintf(stderr,
                  "[ERROR] Failed parsing Clm parameter in line: [%s]; icgem "
                  "file %s (traceback: %s)\n",
                  line, filename.c_str(), __func__);
          ++error;
        }

        ccres = std::from_chars(skip_ws(ccres.ptr), line + sz, Slm);
        if (ccres.ec != std::errc{}) {
          fprintf(stderr,
                  "[ERROR] Failed parsing Slm parameter in line: [%s]; icgem "
                  "file %s (traceback: %s)\n",
                  line, filename.c_str(), __func__);
          ++error;
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
    
    // gfct lines are for tvg field (if L=0=M, no effect ...)
    } else if (starts_with("gfct", line)) {
      // expecting columns: degree, order, Clm, Slm, [...] and time!
      start = line + 4;
      int ll;
      auto ccres = std::from_chars(skip_ws(start), line+sz, ll);
      if (ccres.ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      int mm;
      ccres = std::from_chars(skip_ws(ccres.ptr), line+sz, mm);
      if (ccres.ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing order parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        ++error;
      }

      // only interested in the coefficients, if degree and order are less than
      // max
      if (ll <= l && mm <= m) {

        ccres = std::from_chars(skip_ws(ccres.ptr), line + sz, Clm);
        if (ccres.ec != std::errc{}) {
          fprintf(stderr,
                  "[ERROR] Failed parsing Clm parameter in line: [%s]; icgem "
                  "file %s (traceback: %s)\n",
                  line, filename.c_str(), __func__);
          ++error;
        }

        ccres = std::from_chars(skip_ws(ccres.ptr), line + sz, Slm);
        if (ccres.ec != std::errc{}) {
          fprintf(stderr,
                  "[ERROR] Failed parsing Slm parameter in line: [%s]; icgem "
                  "file %s (traceback: %s)\n",
                  line, filename.c_str(), __func__);
          ++error;
        }

        // do not know exactly at which columns we'll find time, count them
        int cols = count_columns(line);

        // tstart - tstop is found in the last two columns, format is:
        // [yyyymmdd.xxxx]
        dso::datetime<dso::nanoseconds> tstart, tend;
        start = goto_column(line, cols-2);
        if (resolve_date(start, tstart))
          ++error;
        if (resolve_date(start+14. tend))
          ++error;
        
        if (t>= tstart && t < tend) {
          // add drift to harmonic coefficients matrix
          coeffs->C(ll, mm) = Clm;
          ++coeffs_read;
          if (mm == 0) {
            assert(Slm == 0e0);
          } else {
            coeffs->S(ll, mm) = Slm;
          }
        }
      }
    
    // trnd lines are for trend/drift (if L=0=M, no effect ...)
    } else if (starts_with("trnd", line)) {
      // expecting that this 'trnd' field should match the current degree and
      // order of the --already read-- TVG coefficients
      start = line + 4;
      int ll;
      auto ccres = std::from_chars(skip_ws(start), line+sz, ll);
      if (ccres.ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        ++error;
      }

      int mm;
      ccres = std::from_chars(skip_ws(ccres.ptr), line+sz, mm);
      if (ccres.ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing order parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        ++error;
      }

      // there are optional fields, could be missing, skip check
      //if ((ll != tvg_ll) || (mm != tvg_mm)) {
      //  fprintf(stderr,
      //          "[ERROR] Reading line of type \'trnd\' but order/degree do not "
      //          "match with previous TVG coefficients read (%d,%d)!\n",
      //          ll, mm);
      //  fprintf(stderr,
      //          "[ERROR] Current TVG degree and order: %d/%d, icgem file: %s "
      //          "(traceback: %s)\n",
      //          tvg_ll, tvg_mm, filename.c_str(), __func__);
      //  return 1;
      //}
      
      // only interested in the coefficients, if degree and order are less than
      // max
      if (ll <= l && mm <= m) {

        ccres = std::from_chars(skip_ws(ccres.ptr), line + sz, Clm);
        if (ccres.ec != std::errc{}) {
          fprintf(stderr,
                  "[ERROR] Failed parsing Clm parameter in line: [%s]; icgem "
                  "file %s (traceback: %s)\n",
                  line, filename.c_str(), __func__);
          ++error;
        }

        ccres = std::from_chars(skip_ws(ccres.ptr), line + sz, Slm);
        if (ccres.ec != std::errc{}) {
          fprintf(stderr,
                  "[ERROR] Failed parsing Slm parameter in line: [%s]; icgem "
                  "file %s (traceback: %s)\n",
                  line, filename.c_str(), __func__);
          ++error;
        }

        // do not know exactly at which columns we'll find time, count them
        int cols = count_columns(line);

        // tstart - tstop is found in the last two columns, format is:
        // [yyyymmdd.xxxx]
        dso::datetime<dso::nanoseconds> tstart, tend;
        start = goto_column(line, cols-2);
        if (resolve_date(start, tstart))
          ++error;
        if (resolve_date(start+14. tend))
          ++error;
        
        if (t>= tstart && t < tend) {
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
