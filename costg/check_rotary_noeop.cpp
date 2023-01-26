#include "eop.hpp"
#include "geodesy/units.hpp"
#include "orbit_integration.hpp"
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <datetime/dtfund.hpp>
#include <fstream>
#include <vector>
#ifdef USE_SOFA
#include "sofa.h"
#endif

using namespace std::chrono;

dso::TwoPartDate gps2tt(const dso::TwoPartDate &gpst) noexcept {
  constexpr const double offset = (32.184e0 + 19e0) / 86400e0;
  return dso::TwoPartDate(gpst._big, gpst._small + offset).normalized();
}

dso::TwoPartDate gps2tai(const dso::TwoPartDate &gpst) noexcept {
  constexpr const double offset = (19e0) / 86400e0;
  return dso::TwoPartDate(gpst._big, gpst._small + offset).normalized();
}


// gps time [mjd], xp, yp, s' [rad], dUT1, LOD [seconds], X, Y, s [rad]
struct Eop {
  dso::TwoPartDate mjd;
  double xp,yp,sp,dUT1,LOD,X,Y,s;
  dso::EopRecord toEopRecord() const {
    return dso::EopRecord({gps2tt(mjd), dso::rad2sec(xp), dso::rad2sec(yp), dUT1,
                          LOD, dso::rad2sec(X), dso::rad2sec(Y),0e0});
  }
};

// gps time [mjd], xx, xy, xz, yx, yy, yz, zx, zy, zz [-]
struct Rotary {
  dso::TwoPartDate mjd;
  Eigen::Matrix<double,3,3> rot;
};

int map_eopsin(const char *fn, std::vector<Eop> &eops);
int map_rotary(const char *fn, std::vector<Rotary> &rots);
inline const char *skipws(const char *line) noexcept {
  const char *str = line;
  while (*str && *str == ' ')
    ++str;
  return str;
}

int main(int argc, char *argv[]) {
 // check input
  if (argc != 3) {
    fprintf(stderr,
            "USAGE: %s [01earthRotation_interpolatedEOP] [01earthRotation_rotaryMatrix]\n",
            argv[0]);
    return 1;
  }

  // read in input Eops
  std::vector<Eop> eops;
  if (map_eopsin(argv[1], eops)) {
    fprintf(stderr, "ERROR Failed reading file %s\n", argv[1]);
    return 1;
  }

  // read in rotary matrix
  std::vector<Rotary> rots;
  if (map_rotary(argv[2], rots)) {
    fprintf(stderr, "ERROR Failed reading file %s\n", argv[2]);
    return 1;
  }

  Eigen::Matrix<double,3,3> rc2i,rpom;
  double era,xlod;
  std::vector<Rotary>::const_iterator it = rots.cbegin();
  for (auto const &e : eops) {
    // construct TRF->CRF matrix, using parameters of current Eop instance
    if (dso::gcrs2itrs(gps2tai(e.mjd), e.toEopRecord(), e.X, e.Y, e.s, e.sp,
                       rc2i, era, rpom, xlod)) {
      fprintf(stderr,
              "ERROR Failed to compute TRF->CRF transformation matrix\n");
      return 1;
    }

    // actual matrix:
    const auto rc2ti = Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) * rc2i;
    const Eigen::Matrix<double, 3, 3> R = (rpom * rc2ti);

    // report differences
    auto cit = std::find_if(it, rots.cend(), [&](const Rotary &p) {
      return std::abs(
                 p.mjd.diff<dso::DateTimeDifferenceType::FractionalSeconds>(
                     e.mjd)) < 1e-6;
    });

    if (cit != rots.cend()) {
      const Eigen::Matrix<double, 3, 3> RC = cit->rot;
      // differences [-]
      printf("%12.5f %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e "
             "%+.15e\n",
             e.mjd.mjd(), R(0, 0) - RC(0, 0), R(0, 1) - RC(0, 1),
             R(0, 2) - RC(0, 2), R(1, 0) - RC(1, 0),
             R(1, 1) - RC(1, 1), R(1, 2) - RC(1, 2),
             R(2, 0) - RC(2, 0), R(2, 1) - RC(2, 1),
             R(2, 2) - RC(2, 2));
    }
    it = cit;
  }

  return 0;
}

int map_eopsin(const char *fn, std::vector<Eop> &eops) {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR Failed to open file %s\n", fn);
    return 1;
  }

  constexpr const int MAX_LINE = 1024;
  char line[MAX_LINE];

  // two header lines, skip
  for (int i = 0; i < 2; i++)
    fin.getline(line, MAX_LINE);

  double data[9];
  while (fin.getline(line, MAX_LINE)) {
    const char *c = line;
    const int sz = std::strlen(line);
    for (int i = 0; i < 9; i++) {
      auto cres = std::from_chars(skipws(c), line + sz, data[i]);
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
    double it;
    const double ft = std::modf(data[0], &it);
    eops.push_back({dso::TwoPartDate(it, ft), data[1], data[2], data[3],
                    data[4], data[5], data[6], data[7], data[8]});
  }
  return 0;
}

int map_rotary(const char *fn, std::vector<Rotary> &rots) {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR Failed to open file %s\n", fn);
    return 1;
  }

  constexpr const int MAX_LINE = 1024;
  char line[MAX_LINE];

  // five header lines, skip
  for (int i = 0; i < 5; i++)
    fin.getline(line, MAX_LINE);

  double data[10];
  while (fin.getline(line, MAX_LINE)) {
    const char *c = line;
    const int sz = std::strlen(line);
    for (int i = 0; i < 10; i++) {
      auto cres = std::from_chars(skipws(c), line + sz, data[i]);
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
    double it;
    const double ft = std::modf(data[0], &it);
    rots.push_back({
        dso::TwoPartDate(it, ft),
        Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(data + 1)});
  }
  return 0;
}
