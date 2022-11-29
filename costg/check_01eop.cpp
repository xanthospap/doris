#include "eop.hpp"
#include <charconv>
#include <cstdio>
#include <datetime/dtfund.hpp>
#include <fstream>
#include <vector>
#include "geodesy/units.hpp"

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

  printf("%12s %12s %12s %12s %12s %12s %12s\n", "Mjd", "xp('')", "yp('')",
         "dut1 (sec)", "lod (sec)", "X ('')", "Y ('')");

  // ok, now for every MJD in the reference file, check our interpolation
  for (const Eop01Record &eop: refeops) {
    // hold results here
    dso::EopRecord myeop;
    
    // need to transform GPST to TT for our interpolation, also note that 
    // costG uses an order-3 interpolation
    if (eop_lut.interpolate(gps2tt(eop.mjd), myeop, 3)) {
      fprintf(stderr, "ERROR. My interpolation failed!\n");
      return 1;
    }

    // trnaform angular units to arcsec
    const auto reop = eop.to_sec();

    // report results:
    printf("%12.5f %.6f %.6f %.7f %.7f %.6f %.6f\n", eop.mjd, 
      std::abs(reop.xp-myeop.xp),
      std::abs(reop.yp-myeop.yp),
      std::abs(reop.dut1-myeop.dut),
      std::abs(reop.lod-myeop.lod),
      std::abs(reop.X-myeop.dx),
      std::abs(reop.Y-myeop.dy));
  }

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
