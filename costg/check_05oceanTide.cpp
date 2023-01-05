#include "eop.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iau.hpp"
#include "orbit_integration.hpp"
#include "planetpos.hpp"
#include "tides.hpp"
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <datetime/dtcalendar.hpp>
#include <datetime/dtfund.hpp>
#include <fstream>
#include <vector>

using namespace std::chrono;

const int Degree = 80;
const int Order = 80;

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
            "USAGE: %s [eopc04_14_IAU2000.62-now] [oceanTide_FES2014b.potential.iers.txt] [00orbit_icrf.txt] "
            "[04solidEarthTide_icrf.txt]\n",
            argv[0]);
    return 1;
  }

  // parse reference results (gravity acceleration)
  std::vector<Acc> refaccs;
  if (map_input(argv[4], refaccs))
    return 1;

  // parse reference position (ICRF)
  std::vector<Pos> refpos;
  if (map_position(argv[3], refpos))
    return 1;

  // fist date in file as datetime instance
  dso::datetime<dso::nanoseconds> d1(
      dso::modified_julian_day(static_cast<int>(refpos[0].mjd._big)),
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

  // for ICRF-to-ITRF transformations
  Eigen::Matrix<double, 3, 3> rc2i, rpom;
  double era, xlod;

  // read ocean tide file
  std::vector<dso::DoodsonOceanTideConstituent> vdds;
  if (dso::memmap_octide_coefficients(argv[2], vdds, Degree, Order, 3, 1e-11)) {
    fprintf(stderr, "Failed reading inpit file %s!\n", argv[2]);
    return 1;
  }
  //char buf[64];
  //for (const auto &d : vdds) {
  //  printf("%s dC+=%.6e dS+=%.6e dC-=%.6e dS-=%.6e\n", d.doodson_number().str(buf), d.delCp(0,0), d.delSp(0,0), d.delCm(0,0), d.delSm(0,0));
  //}
  //return 2;

  // An ocean tide instance
  dso::OceanTide octide(vdds, 0.3986004415E+15, 0.6378136460E+07, Degree, Order);

  std::vector<Acc>::const_iterator it = refaccs.cbegin();
  for (const auto &pos : refpos) {

    /* ECI to ECEF */
    if (dso::gcrs2itrs(gps2tai(pos.mjd), eop_lut, rc2i, era, rpom, xlod)) {
      fprintf(stderr, "Failed transforming ECI to ECEF\n");
      return 2;
    }

    // transform satellite position vector, icrf-to-itrf
    [[maybe_unused]] const Eigen::Matrix<double, 3, 1> cpos =
        dso::rcel2ter(pos.xyz, rc2i, era, rpom);

    // compute acceleration due to ocean tide (ECEF)
    const auto tt = gps2tt(pos.mjd);
    const auto utc = tt.tt2utc();
    dso::EopRecord eops;
    if (int error; (error = eop_lut.interpolate(utc, eops))) {
      fprintf(stderr, "ERROR. Failed getting EOP values (status: %d)\n", error);
      return error;
    }
    const auto ut1 =
        dso::TwoPartDate(utc._big, utc._small + eops.dut / 86400e0);
    Eigen::Matrix<double, 3, 1> ecef_acc;
    octide.acceleration(tt, ut1, cpos, ecef_acc, Degree, Order);

    // acceleration, ITRF-to-ICRF
    Eigen::Matrix<double, 3, 1> acc = dso::rter2cel(ecef_acc, rc2i, era, rpom);

    // find respective element in reference file
    auto cit = std::find_if(it, refaccs.cend(), [&](const Acc &p) {
      return std::abs(
                 p.mjd.diff<dso::DateTimeDifferenceType::FractionalSeconds>(
                     pos.mjd)) < 1e-3;
    });

    if (cit != refaccs.cend()) {
      // printf("%.12f %+.15f %+.15f %+.15f\n", pos.mjd.mjd(),
      //        acc(0) - cit->a(0), acc(1) - cit->a(1), acc(2) - cit->a(2));
      printf("[MINE] %.12f %+.15e %+.15e %+.15e\n", pos.mjd.mjd(), acc(0),
             acc(1), acc(2));
      printf("[GRPS] %.12f %+.15e %+.15e %+.15e\n", pos.mjd.mjd(), cit->a(0),
             cit->a(1), cit->a(2));
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
