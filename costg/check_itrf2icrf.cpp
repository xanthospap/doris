#include "eop.hpp"
#include "geodesy/units.hpp"
#include "orbit_integration.hpp"
#include <charconv>
#include <cmath>
#include <cstdio>
#include <datetime/dtfund.hpp>
#include <fstream>
#include <vector>
#include <chrono>
using namespace std::chrono;

struct Pos {
  double mjd;
  Eigen::Matrix<double, 3, 1> xyz;
};

int map_position(const char *fn, std::vector<Pos> &rots);
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
double gps2tai(double gpst) noexcept {
  constexpr const double offset = (19e0) / 86400e0;
  return gpst + offset;
}

int main(int argc, char *argv[]) {

  // check input
  if (argc != 4) {
    fprintf(stderr,
            "USAGE: %s [EOP INPUT] [00orbit_itrf.txt] "
            "[00orbit_icrf.txt]\n",
            argv[0]);
    return 1;
  }

  // parse position data
  std::vector<Pos> icrfPos, itrfPos;
  if (map_position(argv[2], itrfPos))
    return 1;
  if (map_position(argv[3], icrfPos))
    return 1;
  
  // fist date in file as datetime instance
  dso::datetime<dso::nanoseconds> d1(
      dso::modified_julian_day(static_cast<int>(itrfPos[0].mjd)),
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
  
  Eigen::Matrix<double, 3, 3> rc2i, rpom;
  double era, xlod;
  
  std::vector<Pos>::const_iterator it = icrfPos.cbegin();
  // for every sattellite position ...
  auto clock_start = high_resolution_clock::now();
  for (const auto &tpos : itrfPos) {

    // Transformation data
    assert(!dso::gcrs2itrs(gps2tai(tpos.mjd), eop_lut, rc2i, era, rpom, xlod));

    // transform itrf 2 icrf
    const Eigen::Matrix<double, 3, 1> cpos =
        dso::rter2cel(tpos.xyz, rc2i, era, rpom);

    // report differences
    auto cit = std::find_if(it, icrfPos.cend(), [&](const Pos &p) {
      return std::abs(p.mjd - tpos.mjd) < 1e-16;
    });
    
    if (cit != icrfPos.cend()) {
      printf("%.12f %+.15f %+.15f %+.15f\n", tpos.mjd, cpos(0) - cit->xyz(0),
             cpos(1) - cit->xyz(1), cpos(2) - cit->xyz(2));
      it = cit;
    }
  }
  
  auto clock_stop = high_resolution_clock::now();

  auto duration = duration_cast<milliseconds>(clock_stop - clock_start);
  std::cout << "Time taken by function: " << duration.count() << " milliseconds\n";

  return 0;
}

// file: 00orbit_i[tc]rf.txt
// Assume columns:
// Time [MJD] x [m] y [m] z [m] x [m/s] y [m/s] z [m/s] acc x [m/s^2] acc y
// [m/s^2] acc z [m/s^2]
int map_position(const char *fn, std::vector<Pos> &poss) {
  constexpr const int MAX_CHARS = 1024;
  char line[MAX_CHARS];
  poss.clear();

  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR. Failed opening file %s\n", fn);
    return 1;
  }

  // first five lines are GROOPS related
  for (int i = 0; i < 5; i++)
    fin.getline(line, MAX_CHARS);

  // read data
  double _data[4];
  while (fin.getline(line, MAX_CHARS)) {
    const char *c = line;
    const int sz = std::strlen(line);
    for (int i = 0; i < 4; i++) {
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
    // remember, COLUMN-WISE order!
    poss.push_back(
        {_data[0], Eigen::Map<Eigen::Matrix<double, 3, 1>>(_data + 1)});
  }

  if (!fin.good() && fin.eof())
    return 0;
  return 3;
}
