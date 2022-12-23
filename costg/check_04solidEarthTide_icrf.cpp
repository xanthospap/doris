#include "eop.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iau.hpp"
#include "orbit_integration.hpp"
#include <charconv>
#include <cmath>
#include <cstdio>
#include <datetime/dtcalendar.hpp>
#include <datetime/dtfund.hpp>
#include <fstream>
#include <vector>
#include <chrono>
#include "planetpos.hpp"
#include "tides.hpp"

using namespace std::chrono;

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
  if (argc != 4) {
    fprintf(stderr,
            "USAGE: %s [eopc04_14_IAU2000.62-now] [00orbit_icrf.txt] "
            "[04solidEarthTide_icrf.txt]\n",
            argv[0]);
    return 1;
  }

  // Load CSPICE/NAIF Kernels
  // -------------------------------------------------------------------------
  int cerror = dso::cspice::load_if_unloaded_spk("data/jpl/de440.bsp");
  cerror += dso::cspice::load_if_unloaded_lsk("data/jpl/naif0012.tls");
  cerror += dso::cspice::load_if_unloaded_pck("data/jpl/gm_de431.tpc");
  // cerror += dso::cspice::load_if_unloaded_pck("data/jpl/earth_fixed.tf");
  // cerror += dso::cspice::load_if_unloaded_pck("data/jpl/earth_latest_high_prec.bpc");
  cerror += dso::cspice::load_if_unloaded_pck("data/jpl/pck00010.tpc");
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

  // get gravitational constant for Sun and Moon, note that these are in 
  // km^3/s^2
  double GMSun, GMMoon;
  if (dso::get_sun_moon_GM("data/jpl/gm_de431.tpc", GMSun, GMMoon)) {
    fprintf(stderr, "Failed getting gravitational constants\n");
    return 1;
  }
  printf("Got Moon/Sun constants from %s GMsun=%.6e GMmoon=%.6e\n",
         "data/jpl/gm_de431.tpc", GMSun, GMMoon);

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

  // A solid earth tide instance (km^3 to m^3)
  dso::SolidEarthTide setide(0.3986004415E+15, 0.6378136460E+07, GMMoon * 1e9,
                             GMSun * 1e9);

  std::vector<Acc>::const_iterator it = refaccs.cbegin();
  for (const auto &pos : refpos) {

    // get Sun and Moon coordinates in ECEF (Cartesian, km)
    double dummy;
    double rm[3], rs[3];
    double et = dso::cspice::jd2et(gps2tt(pos.mjd).mjd() + dso::mjd0_jd);
    
    /* coordinates in km, ECI
    spkezp_c(301, et, "J2000", "NONE", 399, rm, &dummy);
    spkezp_c( 10, et, "J2000", "NONE", 399, rs, &dummy);
    Eigen::Matrix<double,3,1> rmoon(rm[0], rm[1], rm[2]); // [km]
    Eigen::Matrix<double,3,1> rsun(rs[0], rs[1], rs[2]);  // [km]
    rmoon *= 1e3; rsun *= 1e3;
    printf("\tTime is %.20f\n", gps2tt(pos.mjd).mjd());
    printf("\tMoon at %.3f %.3f %.3f\n", rmoon(0), rmoon(1),rmoon(2));
    printf("\tSun at  %.3f %.3f %.3f\n", rsun(0), rsun(1),rsun(2));
    */
    
    /* ECI to ECEF */
    if (dso::gcrs2itrs(gps2tai(pos.mjd), eop_lut, rc2i, era, rpom, xlod)) {
      fprintf(stderr, "Failed transforming ECI to ECEF\n");
      return 2;
    }
    /*
    rmoon = dso::rcel2ter(rmoon, rc2i, era, rpom);
    rsun  = dso::rcel2ter(rsun, rc2i, era, rpom);
    printf("\tMoon at %.3f %.3f %.3f\n", rmoon(0), rmoon(1),rmoon(2));
    printf("\tSun at  %.3f %.3f %.3f\n", rsun(0), rsun(1),rsun(2));
    double sph[3];
    rm[0] = rmoon(0); rm[1]=rmoon(1); rm[2]=rmoon(2);
    recsph_c(rm, &sph[0], &sph[1], &sph[2]);
    printf("\tMoon at %.1f %.15e %.15e\n", sph[0], dso::rad2deg(sph[1]), dso::rad2deg(sph[2]));
    */
 
    /* directly in ECEF from CSPICE */   
    spkezp_c(301, et, "IAU_EARTH", "NONE", 399, rm, &dummy);
    spkezp_c( 10, et, "IAU_EARTH", "NONE", 399, rs, &dummy);
    const Eigen::Matrix<double, 3, 1> rmoon =
        Eigen::Matrix<double, 3, 1>(rm[0], rm[1], rm[2]) * 1e3; // [km]
    const Eigen::Matrix<double, 3, 1> rsun =
        Eigen::Matrix<double, 3, 1>(rs[0], rs[1], rs[2]) * 1e3; // [km]
    //printf("\tMoon at %.3f %.3f %.3f\n", rmoon(0), rmoon(1),rmoon(2));
    //printf("\tSun at  %.3f %.3f %.3f\n", rsun(0), rsun(1),rsun(2));
    //rm[0] = rmoon(0); rm[1]=rmoon(1); rm[2]=rmoon(2);
    //recsph_c(rm, &sph[0], &sph[1], &sph[2]);
    //printf("\tMoon at %.1f %.15e %.15e\n", sph[0], dso::rad2deg(sph[1]), dso::rad2deg(sph[2]));
    
    // transform icrf-to-itrf
    [[maybe_unused]]const Eigen::Matrix<double, 3, 1> cpos =
        dso::rcel2ter(pos.xyz, rc2i, era, rpom);

    // compute acceleration due to solid earth tide (ECEF)
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
    setide.acceleration(tt, ut1, cpos, rmoon, rsun, ecef_acc);

    // acceleration, ITRF-to-ICRF
    Eigen::Matrix<double,3,1> acc = dso::rter2cel(ecef_acc, rc2i, era, rpom);

    // find respective element in reference file
    auto cit = std::find_if(it, refaccs.cend(), [&](const Acc &p) {
      return std::abs(
                 p.mjd.diff<dso::DateTimeDifferenceType::FractionalSeconds>(
                     pos.mjd)) < 1e-3;
    });

    if (cit != refaccs.cend()) {
      //printf("%.12f %+.15f %+.15f %+.15f\n", pos.mjd.mjd(),
      //       acc(0) - cit->a(0), acc(1) - cit->a(1), acc(2) - cit->a(2));
      printf("[MINE] %.12f %+.15f %+.15f %+.15f\n", pos.mjd.mjd(), acc(0),
             acc(1), acc(2));
      printf("[GRPS] %.12f %+.15f %+.15f %+.15f\n", pos.mjd.mjd(), cit->a(0),
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
    accs.push_back(
        {dso::TwoPartDate(it,ft), Eigen::Map<Eigen::Matrix<double, 3, 1>>(_data + 1)});
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
    poss.push_back(
        {dso::TwoPartDate(it,ft), Eigen::Map<Eigen::Matrix<double, 3, 1>>(_data + 1)});
  }

  if (!fin.good() && fin.eof())
    return 0;
  return 3;
}
