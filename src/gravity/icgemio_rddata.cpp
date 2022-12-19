#include "datetime/dtcalendar.hpp"
#include "icgemio.hpp"
#include "iers2010/iersc.hpp"
#include <cassert>
#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <cstring>

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
  while (*str && *str == ' ')
    ++str;
  return str;
}

int count_columns(const char *line) noexcept {
  const char *c = line;
  int cols = 0;
  while (*c) {
    while (*c && *c == ' ')
      ++c;
    if (*c && *c != ' ')
      ++cols;
    while (*c && *c != ' ')
      ++c;
  }
  return cols;
}

const char *goto_column(const char *line, int colnr) noexcept {
  // printf("going to column #%d, of string [%s]\n", colnr, line);
  const char *c = line;
  int cols = 0;
  while (*c) {
    while (*c && *c == ' ')
      ++c;
    if (*c && *c != ' ') {
      if (++cols == colnr + 1) {
        // printf("founf column: [%s]\n", c);
        return c;
      } else {
        // printf("new column: [%s]\n", c);
        ;
      }
    }
    while (*c && *c != ' ')
      ++c;
  }
  return nullptr;
}

int resolve_date(const char *date_str,
                 dso::datetime<dso::nanoseconds> &t) noexcept {
  // expected format: t0[yyyymmdd.xxxx]
  using SecIntType = dso::nanoseconds::underlying_type;
  int year, month, dom, error = 0;
  double fraction, dummy;
  auto rc = std::from_chars(date_str, date_str + 4, year);
  error += (rc.ec != std::errc{});
  rc = std::from_chars(date_str + 4, date_str + 6, month);
  error += (rc.ec != std::errc{});
  rc = std::from_chars(date_str + 6, date_str + 8, dom);
  error += (rc.ec != std::errc{});
  rc = std::from_chars(date_str, date_str + 14, fraction);
  fraction = std::modf(fraction, &dummy);
  const double fnanosec = dso::nanoseconds::sec_factor<double>();
  if (error)
    return error;
  try {
    t = dso::datetime<dso::nanoseconds>(
        dso::year(year), dso::month(month), dso::day_of_month(dom),
        dso::nanoseconds(static_cast<SecIntType>(fnanosec)));
  } catch (std::exception &) {
    fprintf(stderr,
            "[ERROR] Failed to resolve datetime field from string: \"%s\" "
            "(traceback: %s)\n",
            date_str, __func__);
    return 1;
  }
  return 0;
}
} // unnamed namespace

/// @warning coeffs should have already been initialized and allocated with
///          enough memmory to hold the (to-be-) parsed coefficients.
int dso::Icgem::parse_data(int l, int m,
                           const dso::datetime<dso::nanoseconds> &t,
                           dso::HarmonicCoeffs *coeffs) noexcept {

  // clear out coeffs
  coeffs->clear();

  // t in fractional years (needed for TVG terms)
  const double tyears = t.as_fractional_years();

  // pre-compute angular periods for harmonic coefficients (if any are indeed
  // found while parsing)
  double angular_year, sannual, cannual, ssemiannual, csemiannual;
  {
    double iyear, fyear;
    fyear = std::modf(tyears, &iyear);
    if (!fyear && iyear)
      fyear = 1e0;
    const double annual_angle = iers2010::D2PI * fyear;
    sannual = std::sin(annual_angle);
    cannual = std::cos(annual_angle);
    const double semi_annual_angle = 2e0 * iers2010::D2PI * fyear;
    ssemiannual = std::sin(semi_annual_angle);
    csemiannual = std::cos(semi_annual_angle);
    angular_year = fyear;
  }

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

  if (error)
    return error;

  // assign gravity model constants
  coeffs->GM() = earth_gravity_constant;
  coeffs->Re() = radius;
  coeffs->normalized() = this->is_normalized();

  fin.seekg(data_section_pos);

  char line[max_data_line];
  const char *start;
  int ll, mm;
  double Clm, Slm;

  while (fin.getline(line, max_data_line) && !error) {

    const auto sz = std::strlen(line);

    // gfc lines are for static-gravity field (if L=0=M, no effect ...)
    if (starts_with("gfc ", line)) {
      // expecting columns: degree, order, Clm, Slm, [...]; note that it
      // (seldom) happens that the doubles are written in fortran format ...
      start = line + 4;
      auto ccres = std::from_chars(skip_ws(start), line + sz, ll);
      if (ccres.ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        ++error;
      }

      ccres = std::from_chars(skip_ws(ccres.ptr), line + sz, mm);
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
#ifdef DEBUG
        if (ll == 2 && mm == 1)
          printf("Reading coefficient C(2,1)=%.15e\n", Clm);
#endif
        coeffs->C(ll, mm) += Clm;
        if (mm == 0) {
          assert(Slm == 0e0);
        } else {
          coeffs->S(ll, mm) += Slm;
        }
      }

      // gfct lines are for tvg field (if L=0=M, no effect ...)
    } else if (starts_with("gfct", line)) {
      // expecting columns: degree, order, Clm, Slm, [...] and time!
      start = line + 4;
      auto ccres = std::from_chars(skip_ws(start), line + sz, ll);
      if (ccres.ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        return 1;
      }

      ccres = std::from_chars(skip_ws(ccres.ptr), line + sz, mm);
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
        start = goto_column(line, cols - 2);
        if (resolve_date(start, tstart))
          ++error;
        if (resolve_date(start + 14, tend))
          ++error;

        if (t >= tstart && t < tend) {
#ifdef DEBUG
          if (ll == 2 && mm == 1)
            printf("Reading coefficient C(2,1)=%.15e\n", Clm);
#endif
          // add drift to harmonic coefficients matrix
          coeffs->C(ll, mm) += Clm;
          if (mm == 0) {
            assert(Slm == 0e0);
          } else {
            coeffs->S(ll, mm) += Slm;
          }
        }
      }

      // trnd lines are for trend/drift (if L=0=M, no effect ...)
    } else if (starts_with("trnd", line)) {
      start = line + 4;
      auto ccres = std::from_chars(skip_ws(start), line + sz, ll);
      if (ccres.ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        ++error;
      }

      ccres = std::from_chars(skip_ws(ccres.ptr), line + sz, mm);
      if (ccres.ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing order parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        ++error;
      }

      // there are optional fields, could be missing, skip check
      // if ((ll != tvg_ll) || (mm != tvg_mm)) {
      //  fprintf(stderr,
      //          "[ERROR] Reading line of type \'trnd\' but order/degree do not
      //          " "match with previous TVG coefficients read (%d,%d)!\n", ll,
      //          mm);
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
        start = goto_column(line, cols - 2);
        if (resolve_date(start, tstart))
          ++error;
        if (resolve_date(start + 14, tend))
          ++error;

        if (t >= tstart && t < tend) {
          const double dyears = tyears - tstart.as_fractional_years();
          // assign to harmonic coefficients matrix
          coeffs->C(ll, mm) += Clm * dyears;
          if (mm == 0) {
            assert(Slm == 0e0);
          } else {
            coeffs->S(ll, mm) += Slm * dyears;
          }
        }
      }

      // asin/acos lines are for harmonics (if L=0=M, no effect ...)
    } else if (starts_with("asin", line) || starts_with("acos", line)) {
      start = line + 4;
      auto ccres = std::from_chars(skip_ws(start), line + sz, ll);
      if (ccres.ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing degree parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        ++error;
      }

      ccres = std::from_chars(skip_ws(ccres.ptr), line + sz, mm);
      if (ccres.ec != std::errc{}) {
        fprintf(stderr,
                "[ERROR] Failed parsing order parameter in line: [%s]; icgem "
                "file %s (traceback: %s)\n",
                line, filename.c_str(), __func__);
        ++error;
      }

      // there are optional fields, could be missing, skip check
      // if ((ll != tvg_ll) || (mm != tvg_mm)) {
      //  fprintf(stderr,
      //          "[ERROR] Reading line of type \'asin/acos\' but order/degree
      //          do not " "match with previous TVG coefficients read
      //          (%d,%d)!\n", ll, mm);
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
        start = goto_column(line, cols - 3);
        if (resolve_date(start, tstart))
          ++error;
        if (resolve_date(start + 14, tend))
          ++error;
        // find the yearly period, aka annual, semi-annual, etc ...
        double yper;
        ccres = std::from_chars(goto_column(line, cols - 1), line + sz, yper);
        if (ccres.ec != std::errc{})
          ++error;

        if (t >= tstart && t < tend) {
          // angular frequency, 2π/T * delta-years
          if (std::abs(yper - 1e0) < 1e-9) {
            // Annual signal
            coeffs->C(ll, mm) += Clm * cannual;
            if (mm == 0)
              assert(Slm == 0e0);
            else
              coeffs->S(ll, mm) += Slm * sannual;
          } else if (std::abs(yper - 0.5e0) < 1e-9) {
            // Semi-Annual signal
            coeffs->C(ll, mm) += Clm * csemiannual;
            if (mm == 0)
              assert(Slm == 0e0);
            else
              coeffs->S(ll, mm) += Slm * ssemiannual;
          } else {
            // angular frequency
            const double angular_freq = (iers2010::D2PI * angular_year) / yper;
            coeffs->C(ll, mm) += Clm * std::cos(angular_freq);
            if (mm == 0)
              assert(Slm == 0e0);
            else
              coeffs->S(ll, mm) += Slm * std::sin(angular_freq);
          }
        }
      }

    } else {
      printf("Skipping line [%s], don't know what to do with it!\n", line);
    }

    // end resolving line

  } // end reading lines

  // check the status
  if (!fin.good() && fin.eof()) {
    fin.clear();
    return error;
  }

  return -1;
}
