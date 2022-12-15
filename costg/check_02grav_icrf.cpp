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
using namespace std::chrono;

constexpr const int degree = 120;
constexpr const int order = 120;

struct Acc {
  dso::TwoPartDate mjd;
  Eigen::Matrix<double, 3, 1> a;
};

struct Pos {
  dso::TwoPartDate mjd;
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

dso::TwoPartDate gps2tt(const dso::TwoPartDate &gpst) noexcept {
  constexpr const double offset = (32.184e0 + 19e0) / 86400e0;
  return dso::TwoPartDate(gpst._big, gpst._small + offset).normalized();
}
dso::TwoPartDate gps2tai(const dso::TwoPartDate &gpst) noexcept {
  constexpr const double offset = (19e0) / 86400e0;
  return dso::TwoPartDate(gpst._big, gpst._small + offset).normalized();
}

int main(int argc, char *argv[]) {

  // check input
  if (argc != 5) {
    fprintf(stderr,
            "USAGE: %s [gravity field (.gfc)] [EOP INPUT] [00orbit_itrf.txt] "
            "[02gravityfield_icrf.txt]\n",
            argv[0]);
    return 1;
  }

  // parse reference results (gravity acceleration)
  std::vector<Acc> refaccs;
  if (map_input(argv[4], refaccs))
    return 1;

  // parse reference position
  std::vector<Pos> refposs;
  if (map_position(argv[3], refposs))
    return 1;

  // fist date in file as datetime instance
  dso::datetime<dso::nanoseconds> d1(
      dso::modified_julian_day(static_cast<int>(refaccs[0].mjd._big)),
      dso::nanoseconds(0));

  // Parse the input EOP data file to create an EopLookUpTable eop_lut
  dso::EopLookUpTable eop_lut;
  const int ref_mjd = d1.as_mjd();
  const int start = ref_mjd - 5;
  const int end = ref_mjd + 6;
  // parse C04 EOPs and convert time-stamps to TT (not UTC)
  if (parse_iers_C04(argv[2], start, end, eop_lut)) {
    fprintf(stderr, "ERROR. Failed collecting EOP data\n");
    return 1;
  }

  // regularize ERP (DUT and LOD)
  // eop_lut.regularize();

  // handle gravity field and allocate memory
  dso::HarmonicCoeffs harmonics(degree);
  if (dso::parse_gravity_model(argv[1], degree, order, d1, harmonics, false)) {
    fprintf(stderr, "ERROR Reading gfc gravity harmonics file: %s\n", argv[1]);
    return 1;
  }

  // Eigen::Matrix<double, 3, 3> gpartials;
  std::vector<Acc>::const_iterator it = refaccs.cbegin();

  [[maybe_unused]] Eigen::Matrix<double, 3, 1> acc0, acc1, acc2, acc3, acc32,
      acc4;
  [[maybe_unused]] Eigen::Matrix<double, 3, 3> grad;
  [[maybe_unused]] int dummy_it = 0;
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> Mwork(degree + 3,
                                                                degree + 3),
      Wwork(degree + 3, degree + 3);
  Eigen::Matrix<double, 3, 3> rc2i, rpom;
  double era, xlod;

  // for every sattellite position ...
  auto clock_start = high_resolution_clock::now();
  for (const auto &pos : refposs) {
    //if (test::gravacc3(harmonics, pos.xyz, degree, harmonics.Re(),
    //                   harmonics.GM(), acc3, grad))
    //  return 1;
    if (test::gravacc3(harmonics, pos.xyz, degree, harmonics.Re(),
                       harmonics.GM(), acc32, grad, &Wwork, &Mwork))
      return 1;

    // transform acceleration results to ICRF (from ITRF)
    assert(!dso::gcrs2itrs(gps2tai(pos.mjd), eop_lut, rc2i, era, rpom, xlod));

    // find relative reference acceleration result
    auto cit = std::find_if(it, refaccs.cend(), [&](const Acc &a) {
      return std::abs(a.mjd.diff<dso::DateTimeDifferenceType::FractionalDays>(pos.mjd)) < 1e-3;
    });

    // compute differences
    if (cit != refaccs.cend()) {
      //acc3 = dso::rter2cel(acc3, rc2i, era, rpom);
      acc32 = dso::rter2cel(acc32, rc2i, era, rpom);
      //printf("%.12f %+.15f %+.15f %+.15f [31]\n", pos.mjd, acc3(0) - cit->a(0),
      //       acc3(1) - cit->a(1), acc3(2) - cit->a(2));
      printf("%.12f %+.15f %+.15f %+.15f [32]\n", pos.mjd.mjd(), acc32(0) - cit->a(0),
             acc32(1) - cit->a(1), acc32(2) - cit->a(2));
      it = cit;
    }

    // if (pos.mjd-t0 > 0.5) return 0;
    //if (!dummy_it)
    //  return 0;
  }
  auto clock_stop = high_resolution_clock::now();

  auto duration = duration_cast<milliseconds>(clock_stop - clock_start);
  std::cout << "Time taken by function: " << duration.count()
            << " milliseconds\n";

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
    double it;
    const double ft = std::modf(_data[0], &it);
    accs.push_back({dso::TwoPartDate(it, ft),
                    Eigen::Map<Eigen::Matrix<double, 3, 1>>(_data + 1)});
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
    double it;
    const double ft = std::modf(_data[0], &it);
    poss.push_back({dso::TwoPartDate(it, ft),
                    Eigen::Map<Eigen::Matrix<double, 3, 1>>(_data + 1)});
  }

  if (!fin.good() && fin.eof())
    return 0;
  return 3;
}
