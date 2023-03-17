#include "geodesy/units.hpp"
#include "iers2010/cel2ter.hpp"
#include "iers2010/iau.hpp"
#include "orbit_integration.hpp"
#include "sofa.h"
#include <cassert>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <datetime/dtcalendar.hpp>
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

struct Pos {
  dso::TwoPartDate mjd;
  Eigen::Matrix<double, 3, 1> Pxyz;
  Eigen::Matrix<double, 3, 1> Vxyz;
  Eigen::Matrix<double, 3, 1> Axyz;
};

struct EopData {
  dso::TwoPartDate mjd;
  double xp, yp, sp /*[rad]*/, dUT1, LOD /*[seconds]*/, X, Y, s /*[rad]*/;
  dso::EopRecord toEopRecord() const noexcept {
    return dso::EopRecord{gps2tt(mjd), xp, yp, dUT1, LOD, X, Y, 0e0, 0e0};
  }
};

/* WARNING Takes in MJD and returns JD */
void gps2ut1(const dso::TwoPartDate &gpst, double iers_dut1,
             dso::TwoPartDate &ut1) {
  int error = 0;
  // GPS to TAI: TAI - GPS = 19 seconds
  dso::TwoPartDate tai(gpst);
  tai._small += 19e0 / 86400e0;
  tai = tai.jd_split<dso::TwoPartDate::JdSplitMethod::DT>();

  // First, TAI to UTC
  dso::TwoPartDate utc;
  if (iauTaiutc(tai._big, tai._small, &utc._big, &utc._small)) {
    fprintf(stderr, "ERROR call ti iauTaiutc failed\n");
    ++error;
  }

  // UTC to UT1 using IERS ΔUT1
  if (iauUtcut1(utc._big, utc._small, iers_dut1, &ut1._big, &ut1._small)) {
    fprintf(stderr, "ERROR call ti iauUtcut1 failed\n");
    ++error;
  }

  assert(!error);
}

/* WARNING Takes in MJD and returns JD */
void gps2tt(const dso::TwoPartDate &gpst, dso::TwoPartDate &tt) {
  int error = 0;
  // GPS to TAI: TAI - GPS = 19 seconds
  dso::TwoPartDate tai(gpst);
  tai._small += 19e0 / 86400e0;
  tai = tai.jd_split<dso::TwoPartDate::JdSplitMethod::DT>();
  if (iauTaitt(tai._big, tai._small, &tt._big, &tt._small))
    ++error;
  assert(!error);
}

void foo_mine(const dso::TwoPartDate &gpst, const dso::EopRecord &eops,
              double &X, double &Y, double &s, double &sp, double &era,
              Eigen::Matrix<double, 3, 3> &rc2it) {
  const auto tt = gps2tt(gpst);
  iers2010::sofa::xy06(tt, X, Y);
  s = iers2010::sofa::s06(tt, X, Y);
  sp = iers2010::sofa::sp00(tt);
  {
    const auto utc = tt.tt2tai().tai2utc().normalized();
    const auto ut1p =
        dso::TwoPartDate(utc._big, utc._small + eops.dut / dso::sec_per_day);
    const auto ut1 = ut1p.normalized();
    era = iers2010::sofa::era00(ut1);
  }
  const Eigen::Matrix<double, 3, 3> Rc2i = iers2010::sofa::c2ixys(X, Y, s);
  const Eigen::Matrix<double, 3, 3> Rpom =
      iers2010::sofa::pom00(eops.xp, eops.yp, sp);

  Eigen::Matrix<double, 3, 3> gcrf2tirs =
      Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) * Rc2i;
  rc2it = Rpom * gcrf2tirs;
}

void foo_sofa(const dso::TwoPartDate &gpst, const dso::EopRecord &eops,
              double &X, double &Y, double &s, double &sp, double &era,
              double rc2it[3][3]) {
  int error = 0;

  /* TT in JD */
  dso::TwoPartDate tt;
  gps2tt(gpst, tt);
  const double tt1 = tt._big;
  const double tt2 = tt._small;

  /* CIP and CIO, IAU 2006/2000A. */
  iauXy06(tt1, tt2, &X, &Y);
  s = iauS06(tt1, tt2, X, Y);
  sp = iauSp00(tt1, tt2);

  /* Earth rotation angle. */
  dso::TwoPartDate ut1;
  gps2ut1(gpst, eops.dut, ut1); /* UT1 in (quasi-)JD */
  era = iauEra00(ut1._big, ut1._small);

  /* GCRS to CIRS matrix. */
  double rc2i[3][3];
  iauC2ixys(X, Y, s, rc2i);

  /* Form celestial-terrestrial matrix (no polar motion yet). */
  double rc2ti[3][3];
  iauCr(rc2i, rc2ti);
  iauRz(era, rc2ti);

  /* Polar motion matrix (TIRS->ITRS, IERS 2003). */
  double rpom[3][3];
  iauPom00(eops.xp, eops.yp, sp, rpom);

  /* Form celestial-terrestrial matrix (including polar motion). */
  iauRxr(rpom, rc2ti, rc2it);

  /* GCRS-to-ITRS */
  assert(!error);
}

int map_position(const char *fn, std::vector<Pos> &vpos);
int map_eops(const char *fn, std::vector<EopData> &eops);
inline const char *skipws(const char *line) noexcept {
  const char *str = line;
  while (*str && *str == ' ')
    ++str;
  return str;
}

int main(int argc, char *argv[]) {

  // check input
  if (argc != 4) {
    fprintf(
        stderr,
        "USAGE: %s [01earthRotation_interpolatedEOP.txt] [00orbit_itrf.txt] "
        "[00orbit_icrf.txt]\n",
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
  // parse Eop data
  std::vector<EopData> eops;
  if (map_eops(argv[1], eops))
    return 1;

  auto crf_it = icrf.cbegin();
  auto eop_it = eops.cbegin();
  dso::Itrs2Gcrs Rot;

  for (const auto &pt : itrf) {
    /* corresponding GCRF orbit */
    [[maybe_unused]] auto cit =
        std::find_if(crf_it, icrf.cend(), [&](const Pos &p) {
          return std::abs(
                     p.mjd.diff<dso::DateTimeDifferenceType::FractionalDays>(
                         pt.mjd)) < 1e-12;
        });
    /* corresponding (interpolated) EOPs */
    auto eit = std::find_if(eop_it, eops.cend(), [&](const EopData &p) {
      return std::abs(p.mjd.diff<dso::DateTimeDifferenceType::FractionalDays>(
                 pt.mjd)) < 1e-12;
    });

    const auto gpst = pt.mjd;
    double Xm, Ym, sm, spm, eram;
    double Xs, Ys, ss, sps, eras;
    Eigen::Matrix<double,3,3> rc2im;
    double rc2is[3][3];
    foo_mine(gpst, eit->toEopRecord(), Xm, Ym, sm, spm, eram, rc2im);
    foo_sofa(gpst, eit->toEopRecord(), Xs, Ys, ss, sps, eras, rc2is);

    Eigen::Matrix<double,3,1> gcrfm = rc2im.transpose() * pt.Pxyz;
    double gcrfs[3],p[]={pt.Pxyz(0),pt.Pxyz(1),pt.Pxyz(2)};
    iauTrxp(rc2is, p, gcrfs);

    printf(
        "-----------------------------------------------------------------\n");
    printf("[SOFA] %+.18e %+.18e %+.18e %+.18f %+.18f %+.12f %+.12f %+.12f\n", Xs, Ys, ss, sps,
           eras, gcrfs[0], gcrfs[1], gcrfs[2]);
    printf("[MINE] %+.18e %+.18e %+.18e %+.18f %+.18f %+.12f %+.12f %+.12f\n", Xm, Ym, sm, spm,
           eram, gcrfm(0), gcrfm(1), gcrfm(2));
    printf("[COST] %+.18e %+.18e %+.18e %+.18f %+.18f %+.12f %+.12f %+.12f\n", eit->X, eit->Y,
           eit->s, eit->sp, 0e0, cit->Pxyz(0), cit->Pxyz(1), cit->Pxyz(2));
    printf("[DSCG] %+.18e %+.18e %+.18e %+.18f %+.18f %+.12f %+.12f %+.12f\n",
           eit->X - Xs, eit->Y - Ys, eit->s - ss, eit->sp - sps, 0e0,
           cit->Pxyz(0) - gcrfs[0], cit->Pxyz(1) - gcrfs[1],
           cit->Pxyz(2) - gcrfs[2]);
    printf("[DSMN] %+.18e %+.18e %+.18e %+.18f %+.18f %+.12f %+.12f %+.12f\n",
           Xs - Xm, Ys - Ym, ss - sm, sps - spm, eras - eram,
           gcrfs[0] - gcrfm(0), gcrfs[1] - gcrfm(1), gcrfs[2] - gcrfm(2));

    /*
    Rot.prepare_costg(eit->toEopRecord(), eit->s, eit->sp);
    Eigen::Matrix<double,6,1> ref_state_itrf;
    ref_state_itrf.block<3,1>(0,0) = pt.Pxyz;
    ref_state_itrf.block<3,1>(3,0) = pt.Vxyz;
    Eigen::Matrix<double,6,1> mgcrf = Rot.itrf2gcrf(ref_state_itrf);

    printf("%.18e %+.18e %+.18e %+.18e %+.18e %+.18e %+.18e\n", cit->mjd.mjd(),
           cit->Pxyz(0) - mgcrf(0), cit->Pxyz(1) - mgcrf(1),
           cit->Pxyz(2) - mgcrf(2), cit->Vxyz(0) - mgcrf(3),
           cit->Vxyz(1) - mgcrf(4), cit->Vxyz(2) - mgcrf(5));
    */
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
  double _data[10];
  while (fin.getline(line, MAX_CHARS)) {
    const char *c = line;
    const int sz = std::strlen(line);
    for (int i = 0; i < 10; i++) {
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
                    Eigen::Map<Eigen::Matrix<double, 3, 1>>(_data + 1),
                    Eigen::Map<Eigen::Matrix<double, 3, 1>>(_data + 4),
                    Eigen::Map<Eigen::Matrix<double, 3, 1>>(_data + 7)});
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
