#include "eop.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iau.hpp"
#include "orbit_integration.hpp"
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <datetime/dtfund.hpp>
#include <fstream>
#include <vector>

/*
 * Given the 01earthRotation_rotaryMatrix.txt and 00orbit_itrf.txt files,
 * check the results (of transforming itrf to icrf) against the orbit_icrf.txt
 * file
 */

struct Rotary {
  dso::TwoPartDate mjd;
  Eigen::Matrix<double, 3, 3> R;
};
struct Pos {
  dso::TwoPartDate mjd;
  Eigen::Matrix<double, 3, 1> xyz;
};

int map_input(const char* fn, std::vector<Rotary>& rots);
int map_position(const char *fn, std::vector<Pos> &vpos);
inline const char *skipws(const char *line) noexcept {
  const char *str = line;
  while (*str && *str == ' ')
    ++str;
  return str;
}

int main(int argc, char *argv[]) {

  // check input
  if (argc != 4) {
    fprintf(stderr,
            "USAGE: %s [01earthRotation_rotaryMatrix.txt] [00orbit_itrf.txt] "
            "[orbit_icrf.txt]\n",
            argv[0]);
    return 1;
  }

  // parse reference position, itrf
  std::vector<Pos> itrf;
  if (map_position(argv[2], itrf))
    return 1;
  // parse reference position, icrf
  std::vector<Pos> icrf;
  if (map_position(argv[3], icrf))
    return 1;
  // parse Rotation Matrix (itrf-to-icrf)
  std::vector<Rotary> rots;
  if (map_input(argv[1], rots))
    return 1;

  printf("Data: ITRF Pos: %d ICRF Pos: %d Rot. Matrices: %d\n", (int)itrf.size(), (int)icrf.size(), (int)rots.size());

  auto rot_it = rots.cbegin();
  auto crf_it = icrf.cbegin();
  for (const auto &pt : itrf) {
    auto rit = std::find_if(rot_it, rots.cend(), [&](const Rotary &r) {
      return std::abs(r.mjd.diff<dso::DateTimeDifferenceType::FractionalDays>(pt.mjd)) < 1e-6;
    });
    auto cit = std::find_if(crf_it, icrf.cend(), [&](const Pos &p) {
      return std::abs(p.mjd.diff<dso::DateTimeDifferenceType::FractionalDays>(pt.mjd)) < 1e-6;
    });

    // compute differences
    if (rit != rots.cend() && cit != icrf.cend()) {
      const auto crf2 = rit->R.transpose() * pt.xyz;
      printf("%.12f %.5e %.5e %.5e\n", pt.mjd.mjd(), cit->xyz(0)-crf2(0),cit->xyz(1)-crf2(1),cit->xyz(2)-crf2(2));
      rot_it = rit;
      crf_it = cit;
    }
  }

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
    double it;
    const double ft = std::modf(_data[0], &it);
    poss.push_back({dso::TwoPartDate(it, ft),
                    Eigen::Map<Eigen::Matrix<double, 3, 1>>(_data + 1)});
  }

  if (!fin.good() && fin.eof())
    return 0;
  return 3;
}

// Read and parse file 01earthRotation_rotartMatrix.txt
// Assume columns:
// gps time [mjd], xx, xy, xz, yx, yy, yz, zx, zy, zz
int map_input(const char* fn, std::vector<Rotary>& rots)
{
  constexpr const int MAX_CHARS = 1024;
  char line[MAX_CHARS];
  rots.clear();

  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR. Failed opening file %s\n", fn);
    return 1;
  }

  // first five lines are GROOPS related
  for (int i = 0; i < 6; i++)
    fin.getline(line, MAX_CHARS);

  // read data
  double _data[10];
  while (fin.getline(line, MAX_CHARS)) {
    const char* c = line;
    const int sz = std::strlen(line);
    for (int i = 0; i < 10; i++) {
      auto cres = std::from_chars(skipws(c), line + sz, _data[i]);
      if (cres.ec != std::errc {}) {
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
    double it;
    const double ft = std::modf(_data[0], &it);
    rots.push_back(
        { dso::TwoPartDate(it, ft),
            Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(_data + 1) });
  }

  if (!fin.good() && fin.eof())
    return 0;
  return 3;
}
