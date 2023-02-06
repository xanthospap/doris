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
dso::TwoPartDate gps2tt(const dso::TwoPartDate &gpst) noexcept {
  constexpr const double offset = (32.184e0 + 19e0) / 86400e0;
  return dso::TwoPartDate(gpst._big, gpst._small + offset).normalized();
}

struct Rotary {
  dso::TwoPartDate mjd;
  Eigen::Matrix<double, 3, 3> R;
};
struct Pos {
  dso::TwoPartDate mjd;
  Eigen::Matrix<double, 3, 1> xyz;
};
struct RotQ {
  dso::TwoPartDate mjd;
  Eigen::Quaterniond q;
};
struct EopData {
  dso::TwoPartDate mjd;
  double xp, yp, sp /*[rad]*/, dUT1, LOD /*[seconds]*/, X, Y, s /*[rad]*/;
  /* trf to crf */
  Eigen::Matrix<double, 3, 3> R() const {
    // Form the celestial to intermediate-frame-of-date matrix given the CIP
    // X,Y and the CIO locator s.
    const auto rc2i = iers2010::sofa::c2ixys(X, Y, s); // Q matrix
    // call ERA00 to get the ERA rotation angle (need UT1 datetime)
    const dso::TwoPartDate mjd_utc = gps2tt(mjd).tt2utc();
    dso::TwoPartDate ut1 = mjd_utc;
    ut1._small += dUT1 / 86400e0; // add UT1-UTC, interpolated
    const auto era = iers2010::sofa::era00(ut1.normalized());
    // Form the polar motion matrix (W)
    const auto rpom = iers2010::sofa::pom00(xp, yp, sp); // W matrix
    // const auto rc2ti = Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) *
    // rc2i; return (rpom * rc2ti).transpose();
    const auto Rera =
        Eigen::AngleAxisd(-era, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    return (rpom * Rera * rc2i).transpose();
  }
};

int map_input(const char* fn, std::vector<Rotary>& rots);
int map_position(const char *fn, std::vector<Pos> &vpos);
int map_eops(const char *fn, std::vector<EopData> &eops);
int map_roqs(const char *fn, std::vector<RotQ> &rots);
inline const char *skipws(const char *line) noexcept {
  const char *str = line;
  while (*str && *str == ' ')
    ++str;
  return str;
}

int main(int argc, char *argv[]) {

  // check input
  if (argc != 6) {
    fprintf(stderr,
            "USAGE: %s [01earthRotation_rotaryMatrix.txt] [01earthRotation_interpolatedEOP.txt] [01earthRotation_quaternion.txt] [00orbit_itrf.txt] "
            "[00orbit_icrf.txt]\n",
            argv[0]);
    return 1;
  }

  // parse reference position, itrf
  std::vector<Pos> itrf;
  if (map_position(argv[4], itrf))
    return 1;
  // parse reference position, icrf
  std::vector<Pos> icrf;
  if (map_position(argv[5], icrf))
    return 1;
  // parse Rotation Matrix (itrf-to-icrf)
  std::vector<Rotary> rots;
  if (map_input(argv[1], rots))
    return 1;
  // parse Eop data
  std::vector<EopData> eops;
  if (map_eops(argv[2], eops))
    return 1;
  // parse quaternions
  std::vector<RotQ> roqs;
  if (map_roqs(argv[3], roqs))
    return 1;

  printf("Data: ITRF Pos: %d ICRF Pos: %d Rot. Matrices: %d\n", (int)itrf.size(), (int)icrf.size(), (int)rots.size());

  auto rot_it = rots.cbegin();
  auto crf_it = icrf.cbegin();
  auto eop_it = eops.cbegin();
  auto roq_it = roqs.cbegin();

  for (const auto &pt : itrf) {
    auto rit = std::find_if(rot_it, rots.cend(), [&](const Rotary &r) {
      return std::abs(r.mjd.diff<dso::DateTimeDifferenceType::FractionalDays>(pt.mjd)) < 1e-6;
    });
    auto cit = std::find_if(crf_it, icrf.cend(), [&](const Pos &p) {
      return std::abs(p.mjd.diff<dso::DateTimeDifferenceType::FractionalDays>(pt.mjd)) < 1e-6;
    });
    auto eit = std::find_if(eop_it, eops.cend(), [&](const EopData &p) {
      return std::abs(p.mjd.diff<dso::DateTimeDifferenceType::FractionalDays>(pt.mjd)) < 1e-6;
    });
    auto qit = std::find_if(roq_it, roqs.cend(), [&](const RotQ &p) {
      return std::abs(p.mjd.diff<dso::DateTimeDifferenceType::FractionalDays>(pt.mjd)) < 1e-6;
    });

    // compute differences
    if (rit != rots.cend() && cit != icrf.cend() && eit != eops.cend() && qit != roqs.cend()) {
      const auto crf2 = rit->R.transpose() * pt.xyz;
      printf("[ROT] %.12f %.5e %.5e %.5e\n", pt.mjd.mjd(), cit->xyz(0)-crf2(0),cit->xyz(1)-crf2(1),cit->xyz(2)-crf2(2));
      rot_it = rit;
      crf_it = cit;
      const auto crf3 = eit->R() * pt.xyz;
      printf("[EOP] %.12f %.5e %.5e %.5e\n", pt.mjd.mjd(), cit->xyz(0)-crf3(0),cit->xyz(1)-crf3(1),cit->xyz(2)-crf3(2));
      eop_it = eit;
      const auto crf4 = qit->q.conjugate() * pt.xyz;
      printf("[QRT] %.12f %.5e %.5e %.5e\n", pt.mjd.mjd(), cit->xyz(0)-crf4(0),cit->xyz(1)-crf4(1),cit->xyz(2)-crf4(2));
      roq_it = qit;
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

int map_eops(const char *fn, std::vector<EopData> &eops) {
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
