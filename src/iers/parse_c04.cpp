#include "datetime/utcdates.hpp"
#include "eop.hpp"
#include <charconv>
#include <cstdio>
#include <cstring>
#include <fstream>

namespace {
constexpr const std::size_t MAX_LINE_CHARS = 256;

const char *next_num(const char *line) noexcept {
  const char *ptr = line;
  while (*ptr && *ptr == ' ')
    ++ptr;
  return ptr;
}
}// unnamed namespace

/// Example : https://datacenter.iers.org/data/224/eopc04_14_IAU2000.62-now.txt
/// Format of data:
/// FORMAT(3(I4),I7,2(F11.6),2(F12.7),2(F11.6),2(F11.6),2(F11.7),2(F12.6))
/*  Date      MJD      x          y        UT1-UTC       LOD         dX dY x Err
y Err   UT1-UTC Err  LOD Err     dX Err       dY Err "          "           s s
"         "           "          "          s         s            "           "
     (0h UTC)

1962   1   1  37665  -0.012700   0.213000   0.0326338   0.0017230   0.000000
0.000000   0.030000   0.030000  0.0020000  0.0014000    0.004774    0.002000
1962   1   2  37666  -0.015900   0.214100   0.0320547   0.0016690   0.000000
0.000000   0.030000   0.030000  0.0020000  0.0014000    0.004774    0.002000
*/

int dso::parse_iers_C04(const char *c04fn, dso::modified_julian_day start_mjd,
                        dso::modified_julian_day end_mjd,
                        dso::EopLookUpTable &eoptable, bool utc2tt) noexcept {
  // open file
  std::ifstream fin(c04fn);
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Failed opening EOP (C04) file %s (traceback: "
            "%s)\n",
            c04fn, __func__);
    return 1;
  }

  // if needed, resize the EOP table, interval is [start.end)
  int days = end_mjd.as_underlying_type() - start_mjd.as_underlying_type();
  eoptable.clear();
  eoptable.reserve(days);

  char line[MAX_LINE_CHARS];

  dso::EopRecord rec;
  while (fin.getline(line, MAX_LINE_CHARS)) {
    const char *start;
    const char *end;
    long imjd;
    if (*line && (*line != ' ' && *line != '#')) {
      // skip the data (YYYY MM DD) which is 12 chars length
      // get the mjd (note that this is UTC)
      start = next_num(line + 12);
      end = line + std::strlen(line);
      auto fcr = std::from_chars(start, end, imjd);
      if (fcr.ec != std::errc() || fcr.ptr == start) {
        fprintf(stderr,
                "[ERROR] Failed resolving EOP line [%s]. EOP file is %s "
                "(traceback: %s/err=A)\n",
                line, c04fn, __func__);
        return 1;
      }

      const dso::modified_julian_day cmjd(imjd);
      int error = 0;

      // if the date is ok, collect data
      if (cmjd >= start_mjd && cmjd < end_mjd) {
        error = 0;
        // do we need to transform the given UTC date to TT?
        if (utc2tt) {
          rec.mjd = TwoPartDate((double)imjd, 0e0).utc2tt();
        } else {
          rec.mjd = dso::TwoPartDate((double)imjd, 0e0);
        }
        double data[6];
        for (int i = 0; i < 6; i++) {
          start = next_num(fcr.ptr);
          fcr = std::from_chars(start, end, data[i]);
          if (fcr.ec != std::errc() || fcr.ptr == start) {
            ++error;
          }
        }

        if (error) {
          fprintf(stderr,
                  "[ERROR] Failed resolving EOP line [%s]. EOP file is %s "
                  "(traceback: %s/err=B)\n",
                  line, c04fn, __func__);
          return 2;
        }

        // assign (watch the order)
        rec.xp = data[0];
        rec.yp = data[1];
        rec.dut = data[2];
        rec.lod = data[3];
        rec.dx = data[4];
        rec.dy = data[5];

        eoptable.push_back(rec);
      } else if (cmjd >= end_mjd) {
        break;
      }

    } // if (line[0] != ' ') ...
  }   // end reading lines

  // we were supposed to collect:
  return !(days == eoptable.size());
}
