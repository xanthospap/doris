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

// constants from README
const double GM_Sun = 1.32712442076e20;
const double GM_Moon = 0.49028010560e13;

struct Acc {
  dso::TwoPartDate mjd;
  Eigen::Matrix<double, 3, 1> a;
};

struct Pos {
  dso::TwoPartDate mjd;
  Eigen::Matrix<double, 3, 1> xyz;
};

struct RotQ {
  dso::TwoPartDate mjd;
  Eigen::Quaterniond q;
};

int map_input(const char *fn, std::vector<Acc> &rots);
int map_position(const char *fn, std::vector<Pos> &rots);
int map_roqs(const char *fn, std::vector<RotQ> &rots);
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

  int use_input_eop = false;

  // check input
  if (argc < 4) {
    fprintf(stderr,
            "USAGE: %s [eopc04_14_IAU2000.62-now] [00orbit_icrf.txt] "
            "[04solidEarthTide_icrf.txt] <01earthRotation_quaternion.txt>\n",
            argv[0]);
    return 1;
  } else {
    if (argc == 5)
      use_input_eop = true;
  }

  // Load CSPICE/NAIF Kernels
  // -------------------------------------------------------------------------
  int cerror = dso::cspice::load_if_unloaded_spk("data/jpl/de440.bsp");
  cerror += dso::cspice::load_if_unloaded_lsk("data/jpl/naif0012.tls");
  cerror += dso::cspice::load_if_unloaded_pck("data/jpl/gm_de431.tpc");
  // cerror += dso::cspice::load_if_unloaded_pck("data/jpl/pck00010.tpc");
  if (cerror) {
    fprintf(stderr, "Failed to load NAIF kernels!\n");
    return 1;
  }

  // parse reference results (gravity acceleration)
  std::vector<Acc> refaccs;
  if (map_input(argv[3], refaccs))
    return 1;

  // parse reference position (ICRF)
  std::vector<Pos> refpos;
  if (map_position(argv[2], refpos))
    return 1;

  // parse quaternions
  std::vector<RotQ> roqs;
  if (use_input_eop) {
    if (map_roqs(argv[4], roqs))
      return 1;
  }

  // get gravitational constant for Sun and Moon, note that these are in
  // km^3/s^2
  // double GMSun, GMMoon;
  // if (dso::get_sun_moon_GM("data/jpl/gm_de431.tpc", GMSun, GMMoon)) {
  //  fprintf(stderr, "Failed getting gravitational constants\n");
  //  return 1;
  //}
  // printf("Got Moon/Sun constants from %s GMsun=%.6e GMmoon=%.6e\n",
  //       "data/jpl/gm_de431.tpc", GMSun, GMMoon);

  // fist date in file as datetime instance
  dso::datetime<dso::nanoseconds> d1(
      dso::modified_julian_day(static_cast<int>(refpos[0].mjd._big)),
      dso::nanoseconds(0));

  // Parse the input EOP data file to create an EopLookUpTable eop_lut
  dso::EopLookUpTable eop_lut;
  if (!use_input_eop) {
    const int ref_mjd = d1.as_mjd();
    const int start = ref_mjd - 5;
    const int end = ref_mjd + 6;
    // parse C04 EOPs and convert time-stamps to TT (not UTC)
    if (parse_iers_C04(argv[1], start, end, eop_lut)) {
      fprintf(stderr, "ERROR. Failed collecting EOP data\n");
      return 1;
    }
  }

  // for ICRF-to-ITRF transformations
  Eigen::Matrix<double, 3, 3> rc2i, rpom;
  double era, xlod;

  // A solid earth tide instance (km^3 to m^3)
  dso::SolidEarthTide setide(0.3986004415e+15, 0.6378136460e+07, GM_Moon,
                             GM_Sun);

  std::vector<Acc>::const_iterator it = refaccs.cbegin();
  auto roq_it = roqs.cbegin();
  for (const auto &pos : refpos) {

    auto qit = std::find_if(roq_it, roqs.cend(), [&](const RotQ &p) {
      return std::abs(p.mjd.diff<dso::DateTimeDifferenceType::FractionalDays>(
                 pos.mjd)) < 1e-9;
    });
    if (qit != roqs.cend() || !use_input_eop) {

      // get Sun and Moon coordinates in ECEF (Cartesian, km)
      double rm[3], rs[3];
      const double jd = gps2tt(pos.mjd).mjd() + dso::mjd0_jd;
      // position vector of moon, in J2000, [km]
      dso::cspice::j2planet_pos_from(dso::cspice::jd2et(jd), 301, 399, rm);
      dso::cspice::j2planet_pos_from(dso::cspice::jd2et(jd), 10, 399, rs);
      Eigen::Matrix<double, 3, 1> rMon(rm); // [km]
      Eigen::Matrix<double, 3, 1> rSun(rs); // [km]
      // ECI to ECEF, result in [m]
      Eigen::Matrix<double, 3, 1> rmoon, rsun;
      if (!use_input_eop) {
        if (dso::gcrs2itrs(gps2tai(pos.mjd), eop_lut, rc2i, era, rpom, xlod)) {
          fprintf(stderr, "Failed transforming ECI to ECEF\n");
          return 2;
        }
        rmoon = dso::rcel2ter(rMon, rc2i, era, rpom) * 1e3;
        rsun = dso::rcel2ter(rSun, rc2i, era, rpom) * 1e3;
      } else {
        rmoon = (qit->q * rMon) * 1e3;
        rsun = (qit->q * rSun) * 1e3;
      }

      // transform satellite position vector, icrf-to-itrf
      const Eigen::Matrix<double, 3, 1> cpos =
          (use_input_eop) ? (qit->q.conjugate() * pos.xyz)
                          : dso::rcel2ter(pos.xyz, rc2i, era, rpom);

      // compute acceleration due to solid earth tide (ECEF)
      // we need UT1 (and hence dUT1)
      const auto tt = gps2tt(pos.mjd);
      const auto utc = tt.tt2utc();
      dso::EopRecord eops;
      if (int error; (error = eop_lut.interpolate(utc, eops))) {
        fprintf(stderr, "ERROR. Failed getting EOP values (status: %d)\n",
                error);
        return error;
      }
      const auto ut1 =
          dso::TwoPartDate(utc._big, utc._small + eops.dut / 86400e0)
              .normalized();
      Eigen::Matrix<double, 3, 1> ecef_acc;
      setide.acceleration(tt, ut1, cpos, rmoon, rsun, ecef_acc);

      // acceleration, ITRF-to-ICRF
      Eigen::Matrix<double, 3, 1> acc =
          (use_input_eop) ? (qit->q * ecef_acc)
                          : dso::rter2cel(ecef_acc, rc2i, era, rpom);

      // find respective element in reference file
      auto cit = std::find_if(it, refaccs.cend(), [&](const Acc &p) {
        return std::abs(
                   p.mjd.diff<dso::DateTimeDifferenceType::FractionalSeconds>(
                       pos.mjd)) < 1e-3;
      });

      if (cit != refaccs.cend()) {
        printf("[DIFF] %.12f %+.15e %+.15e %+.15e\n", pos.mjd.mjd(),
               acc(0) - cit->a(0), acc(1) - cit->a(1), acc(2) - cit->a(2));
        printf("[MINE] %.12f %+.15e %+.15e %+.15e\n", pos.mjd.mjd(), acc(0),
               acc(1), acc(2));
        printf("[GRPS] %.12f %+.15e %+.15e %+.15e\n", pos.mjd.mjd(), cit->a(0),
               cit->a(1), cit->a(2));
        it = cit;
      }
    }
    roq_it = qit;
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
int map_roqs(const char *fn, std::vector<RotQ> &rots) {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR Failed to open file %s\n", fn);
    return 1;
  }

  constexpr const int MAX_LINE = 1024;
  char line[MAX_LINE];

  // two header lines, skip
  for (int i = 0; i < 5; i++)
    fin.getline(line, MAX_LINE);

  double data[5];
  while (fin.getline(line, MAX_LINE)) {
    const char *c = line;
    const int sz = std::strlen(line);
    for (int i = 0; i < 5; i++) {
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
    rots.push_back({dso::TwoPartDate(it, ft), Eigen::Quaternion(data[1], data[2], data[3],
                    data[4])});
  }
  return 0;
}
