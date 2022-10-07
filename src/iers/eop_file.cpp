#include "eop.hpp"
#include <cstdio>
#include <cstring>
#include <fstream>

constexpr const std::size_t MAX_LINE_CHARS = 512;

dso::EopFile::EopFile(const char *fn) {
  assert(std::strlen(fn) < 256);
  std::strcpy(filename, fn);
}

dso::EopFile::EopFile(dso::EopFile &&other) noexcept {
  std::strcpy(filename, other.filename);
}

dso::EopFile &dso::EopFile::operator=(dso::EopFile &&other) noexcept {
  std::strcpy(filename, other.filename);
  return *this;
}

int dso::EopFile::parse(dso::modified_julian_day start_mjd,
                        dso::modified_julian_day end_mjd,
                        dso::EopLookUpTable &etable) noexcept {
  // open file
  std::ifstream fin(filename);
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Failed opening Eop file %s (traceback: "
            "%s)\n",
            filename, __func__);
    return -1;
  }

  // if needed, resize the EOP table
  int days = end_mjd.as_underlying_type() - start_mjd.as_underlying_type();
  etable.resize(days);

  char line[MAX_LINE_CHARS];
  int sz = 0;
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
        etable.mjd[sz] = static_cast<double>(mjd);

        start = end;
        etable.xpa[sz] = std::strtod(start, &end);
        error += (start == end);

        start = end;
        etable.ypa[sz] = std::strtod(start, &end);
        error += (start == end);

        start = end;
        etable.ut1a[sz] = std::strtod(start, &end);
        error += (start == end);
        
        start = end;
        etable.dxa[sz] = std::strtod(start, &end);
        error += (start == end);

        start = end;
        etable.dya[sz] = std::strtod(start, &end);
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
  return !(days == sz);
}
