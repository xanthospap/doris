#include "eop.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iau.hpp"
#include "orbit_integration.hpp"
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <datetime/dtcalendar.hpp>
#include <datetime/dtfund.hpp>
#include <fstream>
#include <vector>

struct Acc {
  dso::TwoPartDate mjd;
  Eigen::Matrix<double, 3, 1> a;
};

struct Pos {
  dso::TwoPartDate mjd;
  Eigen::Matrix<double, 3, 1> xyz;
};

int map_position(const char *fn, std::vector<Pos> &vpos);
int map_input(const char *fn, std::vector<Acc> &vacc);
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
            "USAGE: %s [00orbit_icrf.txt] "
            "[03directTideMoon_icrf.txt] [spk kernel, de421] [lsk kernel, naif0012.tls]\n",
            argv[0]);
    return 1;
  }
  
  // parse reference position (icrf)
  std::vector<Pos> icrf;
  if (map_position(argv[1], icrf))
    return 1;

  // parse reference results (gravity acceleration)
  std::vector<Acc> accs;
  if (map_input(argv[2], accs))
    return 1;

  // load NAIF kernels
  dso::cspice::load_if_unloaded_spk(argv[3]);
  dso::cspice::load_if_unloaded_lsk(argv[4]);
  
  Eigen::Matrix<double, 3, 1> moon_acc;
  Eigen::Matrix<double, 3, 3> partials;
  auto acc_it = accs.cbegin();
  for (const auto &crf : icrf) {
    double rmon[3];
    const double jd = gps2tt(crf.mjd).mjd() + dso::mjd0_jd;

    // position vector of moon, in J2000, [km]
    dso::cspice::j2planet_pos_from(dso::cspice::jd2et(jd), 301, 399, rmon);

    // Moon-induced acceleration [m/sec^2]
    Eigen::Matrix<double, 3, 1> rMon(rmon); // [km]
    moon_acc =
        dso::point_mass_accel(0.49028010560e13, crf.xyz, rMon * 1e3, partials);
    
    auto ait = std::find_if(acc_it, accs.cend(), [&](const Acc &p) {
      return std::abs(p.mjd.diff<dso::DateTimeDifferenceType::FractionalDays>(crf.mjd)) < 1e-6;
    });

    if (ait != accs.cend()) {
      printf("%.12f %.5e %.5e %.5e\n", crf.mjd.mjd(), ait->a(0)-moon_acc(0),ait->a(1)-moon_acc(1),ait->a(2)-moon_acc(2));
      // printf("%.12f %.15e %.15e %.15e\n", crf.mjd.mjd(), moon_acc(0),moon_acc(1),moon_acc(2));
      acc_it = ait;
    }
  }

  return 0;
}

// file: 03directTideMoon_icrf.txt
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
  for (int i = 0; i < 6; i++)
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
  for (int i = 0; i < 6; i++)
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
