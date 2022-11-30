#include "eop.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iau.hpp"
#include <charconv>
#include <cmath>
#include <cstdio>
#include <datetime/dtfund.hpp>
#include <fstream>
#include <vector>

class RunningStats {
  // see https://www.johndcook.com/blog/standard_deviation/
private:
  int n;
  double m_old, m_new, s_old, s_new;

public:
  void update(double val) noexcept {
    ++n;
    if (n == 1) {
      m_old = m_new = val;
      m_old = 0e0;
    } else {
      m_new = m_old + (val - m_old) / n;
      s_new = s_old + (val - m_old) * (val - m_new);
      m_old = m_new;
      s_old = s_new;
    }
    return;
  }

  double mean() const noexcept { return (n > 0) ? m_new : 0e0; }
  double variance() const noexcept { return ((n > 1) ? m_new / (n - 1) : 0e0); }
  double stddev() const noexcept { return std::sqrt(variance()); }
};

struct Eop01Record {
  double mjd,         // [mjd]
      xp, yp, sprime, // [rad]
      dut1, lod,      // [sec]
      X, Y, s;        // [rad]

  Eop01Record to_sec() const noexcept {
    Eop01Record n(*this);
    n.xp = dso::rad2sec(xp);
    n.yp = dso::rad2sec(yp);
    n.sprime = dso::rad2sec(sprime);
    n.X = dso::rad2sec(X);
    n.Y = dso::rad2sec(Y);
    n.s = dso::rad2sec(s);
    return n;
  }
};

int map_input(const char *fn, std::vector<Eop01Record> &eops);
inline const char *skipws(const char *line) noexcept {
  const char *str = line;
  while (*str && *str == ' ')
    ++str;
  return str;
}

double gps2tt(double gpst) noexcept {
  constexpr const double offset = (32.184e0 + 19e0) / 86400e0;
  return gpst + offset;
}

int main(int argc, char *argv[]) {

  // check input
  if (argc != 3) {
    fprintf(stderr,
            "USAGE: %s [EOP INPUT] [01earthRotation_interpolatedEOP.txt]\n",
            argv[0]);
    return 1;
  }

  // parse reference results
  std::vector<Eop01Record> refeops;
  if (map_input(argv[2], refeops))
    return 1;

  // fist date in file as datetime instance
  dso::datetime<dso::nanoseconds> d1(
      dso::modified_julian_day(static_cast<int>(refeops[0].mjd)),
      dso::nanoseconds(0));

  // Parse the input EOP data file to create an EopLookUpTable eop_lut
  dso::EopLookUpTable eop_lut;
  const int ref_mjd = d1.as_mjd();
  const int start = ref_mjd - 5;
  const int end = ref_mjd + 6;
  // parse C04 EOPs and convert time-stamps to TT (not UTC)
  if (parse_iers_C04(argv[1], start, end, eop_lut)) {
    fprintf(stderr, "ERROR. Failed collecting EOP data\n");
    return 1;
  }

  // regularize ERP (DUT and LOD)
  eop_lut.regularize();

  // Mean and std. deviation for all parameters
  RunningStats rsxp, rsyp, rsdut1, rslod, rsX, rsY, rsS, rsSprime;

  printf("#%12s %9s %9s %10s %10s %9s %9s %9s %9s\n", "Mjd", "xp('')", "yp('')",
         "dut1 (sec)", "lod (sec)", "X ('')", "Y ('')", "CIO ('')", "TIO ('')");

  // ok, now for every MJD in the reference file, check our interpolation
  for (const Eop01Record &eop : refeops) {
    // hold results here
    dso::EopRecord myeop;

    // need to transform GPST to TT for our interpolation, also note that
    // costG uses an order-3 interpolation
    if (eop_lut.interpolate(gps2tt(eop.mjd), myeop, 3)) {
      fprintf(stderr, "ERROR. My interpolation failed!\n");
      return 1;
    }

    // compute X,Y from series, IAU2006/2000A
    double X, Y;
    iers2010::sofa::xy06(dso::mjd0_jd, gps2tt(eop.mjd), X, Y);
    // compute CIO locator, s [radians]
    const double s = iers2010::sofa::s06(dso::mjd0_jd, gps2tt(eop.mjd), X, Y);
    // compute TIO locator, s' [radians]
    const double sp = iers2010::sofa::sp00(dso::mjd0_jd, gps2tt(eop.mjd));

    // add corrections (from EOP interpolation) to X,Y
    X += dso::sec2rad(myeop.dx);
    Y += dso::sec2rad(myeop.dy);

    // transform angular units to arcsec
    const auto reop = eop.to_sec();

    // report results:
    printf(
  #ifdef VISUAL
        "%12.5f %+.6f %+.6f %+.7f %+.7f %+.6f %+.6f %+.6f %+.6f\n",
  #else
        "%12.5f %+.12e %+.12e %+.12e %+.12e %+.12e %+.12e %+.12e %+.12e\n",
#endif
        eop.mjd, std::abs(reop.xp - myeop.xp), std::abs(reop.yp - myeop.yp),
        std::abs(reop.dut1 - myeop.dut), std::abs(reop.lod - myeop.lod),
        std::abs(dso::rad2sec(eop.X - X)), std::abs(dso::rad2sec(eop.Y - Y)),
        std::abs(dso::rad2sec(eop.s - s)),
        std::abs(dso::rad2sec(eop.sprime - sp)));

    // update statistics
    rsxp.update(reop.xp - myeop.xp);      // arcsec
    rsyp.update(reop.yp - myeop.yp);      // arcsec
    rsdut1.update(reop.dut1 - myeop.dut); // sec
    rslod.update(reop.lod - myeop.lod);   // sec
    rsX.update(eop.X - X);                // rad
    rsY.update(eop.Y - Y);                // rad
    rsS.update(eop.s - s);                // rad
    rsSprime.update(eop.sprime - sp);     // rad
  }

  // print statistics
  printf("#%12s %9s %9s %10s %10s %9s %9s %9s %9s\n", "Mjd", "xp('')", "yp('')",
         "dut1 (sec)", "lod (sec)", "X ('')", "Y ('')", "CIO ('')", "TIO ('')");
  printf("#%12s %+.6f %+.6f %+.7f %+.7f %+.6f %+.6f %+.6f %+.6f\n", " ",
         rsxp.mean(), rsyp.mean(), rsdut1.mean(), rslod.mean(),
         dso::rad2sec(rsX.mean()), dso::rad2sec(rsY.mean()),
         dso::rad2sec(rsS.mean()), dso::rad2sec(rsSprime.mean()));
  printf("#%12s %.6f %.6f %.7f %.7f %.6f %.6f %.6f %.6f\n", " ", rsxp.stddev(),
         rsyp.stddev(), rsdut1.stddev(), rslod.stddev(),
         dso::rad2sec(rsX.stddev()), dso::rad2sec(rsY.stddev()),
         dso::rad2sec(rsS.stddev()), dso::rad2sec(rsSprime.stddev()));

  return 0;
}

// Read and parse file 01earthRotation_interpolatedEOP.txt
// Assume columns:
// gps time [mjd], xp, yp, s' [rad], dUT1, LOD [seconds], X, Y, s [rad]
int map_input(const char *fn, std::vector<Eop01Record> &eops) {
  constexpr const int MAX_CHARS = 1024;
  char line[MAX_CHARS];
  eops.clear();

  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR. Failed opening file %s\n", fn);
    return 1;
  }

  // first two lines are GROOPS related
  fin.getline(line, MAX_CHARS);
  fin.getline(line, MAX_CHARS);

  // read data
  double _data[9];
  while (fin.getline(line, MAX_CHARS)) {
    const char *c = line;
    const int sz = std::strlen(line);
    for (int i = 0; i < 9; i++) {
      auto cres = std::from_chars(skipws(c), line + sz, _data[i]);
      if (cres.ec != std::errc{}) {
        fprintf(stderr, "ERROR Failed resolving line %s\n", line);
        if (!std::strlen(line)) {
          fprintf(stderr, "Line is actually empty, so pretending this never "
                          "happened ...\n");
        } else {
          return 2;
        }
      }
      c = cres.ptr;
    }
    eops.push_back({_data[0], _data[1], _data[2], _data[3], _data[4], _data[5],
                    _data[6], _data[7], _data[8]});
  }

  printf("Number of records collected to compare: %d\n", (int)eops.size());

  if (!fin.good() && fin.eof())
    return 0;
  return 3;
}
