#include "satellites.hpp"
#include <charconv>
#include <cstdio>
#include <fstream>

int dso::get_satellite_corrections(const char *j3mass,
                                   const dso::datetime<dso::nanoseconds> &t,
                                   double *dmdxyz) noexcept {
  std::ifstream fin(j3mass);
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Failed to open satellite information file: %s (traceback: "
            "%s)\n",
            j3mass, __func__);
    return 1;
  }

  const double jd = t.mjd().to_julian_day();
  const double sec = t.sec().to_fractional_seconds();

  for (int i = 0; i < 4; i++)
    dmdxyz[i] = 0e0;

  char line[256];
  int error = 0;
  while (fin.getline(line, 256) && !error) {
    // check for comment lines ...
    if (line[0] != 'C' && line[0] != '/') {
      const char *c = line;
      char *last = line + 255;

      // remove whitespace chars
      while (*c++ == ' ')
        ;

      // first resolve datetime from Julian date
      double cjd, csec;
      auto pec = std::from_chars(c, last, cjd);
      c = pec.ptr;
      if (*c++ != ' ' || pec.ec != std::errc{}) {
        ++error;
      }
      while (*c++ == ' ')
        ;
      pec = std::from_chars(c, last, csec);
      c = pec.ptr;
      if (*c++ != ' ' || pec.ec != std::errc{}) {
        ++error;
      }

      // reached a date later than the requested, end here
      if (error || (cjd >= jd && csec >= sec))
        break;

      // prior date, update corrections
      for (int i = 0; i < 4; i++) {
        while (*c++ == ' ')
          ;
        pec = std::from_chars(c, last, dmdxyz[i]);
        c = pec.ptr + 1;
        if (/**c++ != ' ' || */ pec.ec != std::errc{}) {
          ++error;
        }
      }
    }
  }

  if (error) {
    fprintf(
        stderr,
        "[ERROR] Failed resolving line: \"%s\" in satellite ninformation file: "
        "%s (traceback: "
        "%s)\n",
        line, j3mass, __func__);
    return error;
  }

  return 0;
}
