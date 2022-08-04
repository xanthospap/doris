#include "iers_bulletin.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <datetime/dtfund.hpp>

constexpr const std::size_t MAX_LINE_CHARS = 512;

// retrun a pointer to the first non-whitespace char
inline char *nws(char *line) noexcept {
  char *c = line;
  while (*c && *c == ' ')
    ++c;
  return c;
}

dso::EopFile::EopFile(const char *fn) {
  assert(std::strlen(fn) < 256);
  std::strcpy(filename, fn);
  //if (!stream.is_open()) {
  //  fprintf(stderr,
  //          "[ERROR] Failed opening Eop file %s (traceback: %s)\n",
  //          filename, __func__);
  //  throw std::runtime_error("Failed opening file");
  //}
}

dso::EopFile::EopFile(dso::EopFile &&other) noexcept {
  std::strcpy(filename, other.filename);
}

dso::EopFile &
dso::EopFile::operator=(dso::EopFile &&other) noexcept {
  std::strcpy(filename, other.filename);
  return *this;
}

int dso::EopFile::parse(dso::modified_julian_day start_mjd,
                        dso::modified_julian_day end_mjd, double *mjda, double *xpa,
                        double *ypa, double *ut1a, int &sz) noexcept {
  // no data yet ...
  sz = 0;

  // open file
  std::ifstream fin(filename);
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Failed opening Eop file %s (traceback: "
            "%s)\n",
            filename, __func__);
    return -1;
  }

  char line[MAX_LINE_CHARS];

  // skip any line starting with '#'
  while (fin.getline(line, MAX_LINE_CHARS)) {
    const char *start;
    char *end;
    if (*line != '#') {
      // DATE     MJD       x       y      UT1-UTC      dX     dY     x err    y
      // err   UT1 err  dX err  dY err (0 h UTC)         mas     mas       ms
      // mas mas     mas      mas      ms     mas     mas

      // skip the data (YYYY MM DD) which is 12 chars length
      // get the mjd
      long mjd = std::strtol(line + 12, &end, 10);
      if (mjd == 0 || end == line + 12) {
        fprintf(stderr,
                "[ERROR] Failed resolving EOP line [%s]. Eop file is %s "
                "(traceback: %s)\n",
                line, filename, __func__);
        return -1;
      }

      const dso::modified_julian_day cmjd(mjd);
      int error = 0;

      // if the date is ok, collect data
      if (cmjd >= start_mjd && cmjd < end_mjd) {
        error = 0;
        mjda[sz] = static_cast<double>(mjd);

        start = end;
        xpa[sz] = std::strtod(start, &end);
        error += (start == end);

        start = end;
        ypa[sz] = std::strtod(start, &end);
        error += (start == end);

        start = end;
        ut1a[sz] = std::strtod(start, &end);
        error += (start == end);

        ++sz;
      } else if (cmjd >= end_mjd) {
        break;
      }

      if (error) {
        fprintf(stderr,
                "[ERROR] Failed resolving EOP line [%s]. Eop file is %s "
                "(traceback: %s)\n",
                line, filename, __func__);
        return -1;
      }
    }
  }

  // we were supposed to collect:
  int sz_check = end_mjd.as_underlying_type() - start_mjd.as_underlying_type();
  return !(sz_check == sz);
}
