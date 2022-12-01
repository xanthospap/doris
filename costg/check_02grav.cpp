#include "eop.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iau.hpp"
#include "orbit_integration.hpp"
#include <charconv>
#include <cmath>
#include <cstdio>
#include <datetime/dtfund.hpp>
#include <fstream>
#include <vector>

constexpr const int degree = 90;
constexpr const int order = 90;

struct Acc {
  double mjd;
  Eigen::Matrix<double, 3, 1> a;
};

struct Pos {
  double mjd;
  Eigen::Matrix<double, 3, 1> xyz;
};

int map_input(const char *fn, std::vector<Acc> &rots);
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
            "USAGE: %s [gravity field (.gfc)] [00orbit_itrf.txt] "
            "[02gravityfield_itcrf.txt]\n",
            argv[0]);
    return 1;
  }

  // parse reference results (gravity acceleration)
  std::vector<Acc> refaccs;
  if (map_input(argv[3], refaccs))
    return 1;

  // parse reference position
  std::vector<Pos> refposs;
  if (map_position(argv[2], refposs))
    return 1;

  // fist date in file as datetime instance
  dso::datetime<dso::nanoseconds> d1(
      dso::modified_julian_day(static_cast<int>(refaccs[0].mjd)),
      dso::nanoseconds(0));

  // handle gravity field and allocate memory
  dso::HarmonicCoeffs harmonics(degree);
  if (dso::parse_gravity_model(argv[1], degree, order, d1, harmonics, true)) {
    fprintf(stderr, "ERROR Reading gfc gravity harmonics file: %s\n", argv[1]);
    return 1;
  }
  dso::Mat2D<dso::MatrixStorageType::Trapezoid> V(degree + 3, order + 3),
      W(degree + 3, order + 3);

  Eigen::Matrix<double, 3, 3> gpartials;
  std::vector<Acc>::const_iterator it = refaccs.cbegin();

  // for every sattellite position ...
  for (const auto &pos : refposs) {
    // compute gravity acceleration at this point
    Eigen::Matrix<double, 3, 1> gacc = dso::grav_potential_accel(
        pos.xyz, degree, order, V, W, harmonics, gpartials);
    
    // find relative reference acceleration result
    auto cit = std::find_if(it, refaccs.cend(), [&](const Acc &a) {
      return std::abs(a.mjd - pos.mjd) < 1e-16;
    });
    
    // compute differences
    if (cit != refaccs.cend()) {
      //printf("comparing %.6f %.6f %.6f %.6f\n", pos.mjd, gacc(0),gacc(1),gacc(2));
      //printf("to        %.6f %.6f %.6f %.6f\n", cit->mjd, cit->a(0), cit->a(1), cit->a(2));
      printf("%.12f %.15e %.15e %.15e\n", pos.mjd, gacc(0) - cit->a(0),
             gacc(1) - cit->a(1), gacc(2) - cit->a(2));
      it = cit;
    }
  }

  return 0;
}

// file: 02gravityfield_i[tc]rf.txt
// Assume columns:
// Time [MJD] x [m/s^2] y [m/s^2] z [m/s^2]
int map_input(const char *fn, std::vector<Acc> &accs) {
  constexpr const int MAX_CHARS = 1024;
  char line[MAX_CHARS];
  accs.clear();

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
    accs.push_back(
        {_data[0], Eigen::Map<Eigen::Matrix<double, 3, 1>>(_data + 1)});
  }

  if (!fin.good() && fin.eof())
    return 0;
  return 3;
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
