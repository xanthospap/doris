#include "datetime/utcdates.hpp"
#include "eigen3/Eigen/Geometry"
#include "eop.hpp"
#include "iers2010/iau.hpp"
#include "iers2010/iers2010.hpp"
#include "orbit_integration.hpp"
#include <datetime/dtcalendar.hpp>
#include <iers2010/iersc.hpp>

namespace {
const double TAI2TTFDAYS = 32.184e0 / 86400e0;
}

inline double OmegaEarth(double xlod) noexcept {
  // return iers2010::OmegaEarth * (1e0 - xlod / 86400e0);
  // xlod in [seconds/day]
  // see https://www.iers.org/IERS/EN/Science/EarthRotation/UT1LOD.html
  // tranform LOD to milliseconds (from seconds)
  const double LOD = xlod * 1e3;
  return (72921151.467064e0 - 0.843994809e0 * LOD) * 1e-12; // [rad/sec]
}

Eigen::Matrix<double, 3, 1>
dso::rcel2ter(const Eigen::Matrix<double, 3, 1> r,
              const Eigen::Matrix<double, 3, 3> &rc2i, const double era,
              const Eigen::Matrix<double, 3, 3> &rpom) noexcept {
  const auto rc2ti = Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) * rc2i;
  return (rpom * rc2ti) * r;
}

Eigen::Matrix<double, 6, 1>
dso::ycel2ter(const Eigen::Matrix<double, 6, 1> y,
              const Eigen::Matrix<double, 3, 3> &rc2i, const double era,
              double xlod, const Eigen::Matrix<double, 3, 3> &rpom) noexcept {
  Eigen::Matrix<double, 6, 1> x;
  x.block<3, 1>(0, 0) = dso::rcel2ter(y.block<3, 1>(0, 0), rc2i, era, rpom);
  x.block<3, 1>(3, 0) = dso::rcel2ter(y.block<3, 1>(3, 0), rc2i, era, rpom);

  const double data[] = {0e0, -1e0, 0e0, 1e0, 0e0, 0e0, 0e0, 0e0, 0e0};

  const auto rc2ti = (OmegaEarth(xlod) * Eigen::Matrix<double, 3, 3>(data) *
                      Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ())) *
                     rc2i;
  x.block<3, 1>(3, 0) += (rpom * rc2ti) * y.block<3, 1>(0, 0);
  return x;
}

Eigen::Matrix<double, 3, 1>
dso::rter2cel(const Eigen::Matrix<double, 3, 1> r,
              const Eigen::Matrix<double, 3, 3> &rc2i, const double era,
              const Eigen::Matrix<double, 3, 3> &rpom) noexcept {
  const auto rc2ti = Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) * rc2i;
  return (rpom * rc2ti).transpose() * r;
}

Eigen::Matrix<double, 6, 1>
dso::yter2cel(const Eigen::Matrix<double, 6, 1> y,
              const Eigen::Matrix<double, 3, 3> &rc2i, const double era,
              double xlod, const Eigen::Matrix<double, 3, 3> &rpom) noexcept {
  Eigen::Matrix<double, 6, 1> x;
  x.block<3, 1>(0, 0) = dso::rter2cel(y.block<3, 1>(0, 0), rc2i, era, rpom);
  x.block<3, 1>(3, 0) = dso::rter2cel(y.block<3, 1>(3, 0), rc2i, era, rpom);

  const double data[] = {0e0, -1e0, 0e0, 1e0, 0e0, 0e0, 0e0, 0e0, 0e0};

  const auto rc2ti = (OmegaEarth(xlod) * Eigen::Matrix<double, 3, 3>(data) *
                      Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ())) *
                     rc2i;
  x.block<3, 1>(3, 0) += (rpom * rc2ti).transpose() * y.block<3, 1>(0, 0);
  return x;
}

// IAU 2006/2000A, CIO based, using X,Y series
int dso::gcrs2itrs(const dso::TwoPartDate &mjd_tai,
                   const dso::EopLookUpTable &eop_table,
                   Eigen::Matrix<double, 3, 3> &rc2i, double &era,
                   Eigen::Matrix<double, 3, 3> &rpom, double &xlod) noexcept {

  // TAI to TT
  // const dso::TwoPartDate mjd_tt(mjd_tai._big, mjd_tai._small + TAI2TTFDAYS);
  const dso::TwoPartDate mjd_tt = mjd_tai.tai2tt();

  // TAI to UTC
  const dso::TwoPartDate mjd_utc = mjd_tai.tai2utc();

  // interpolate/correct EOP values using UTC
  dso::EopRecord eops;
  if (int error; (error = eop_table.interpolate(mjd_utc, eops))) {
    fprintf(stderr, "ERROR. Failed getting EOP values (status: %d)\n", error);
    return error;
  }

  // assign interpolated LOD value, [sec/day]
  xlod = eops.lod;

  // split date as in SOFA Date&Time idiom
  // auto sofajd = mjd_tt.jd_sofa();

  // X,Y coordinates of celestial intermediate pole from series based
  // on IAU 2006 precession and IAU 2000A nutation.
  double X, Y;
  iers2010::sofa::xy06(/*sofajd._big, sofajd._small*/mjd_tt, X, Y);

  // The CIO locator s, positioning the Celestial Intermediate Origin on
  // the equator of the Celestial Intermediate Pole, given the CIP's X,Y
  // coordinates. Compatible with IAU 2006/2000A precession-nutation.
  const double s = iers2010::sofa::s06(/*sofajd._big, sofajd._small*/mjd_tt, X, Y);

  // Add CIP corrections (arcsec to radians)
  X += dso::sec2rad(eops.dx);
  Y += dso::sec2rad(eops.dy);

  // Form the celestial to intermediate-frame-of-date matrix given the CIP
  // X,Y and the CIO locator s.
  // [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]  = RC2T * [CRS]
  rc2i = iers2010::sofa::c2ixys(X, Y, s);

  // call ERA00 to get the ERA rotation angle (need UT1 datetime)
  dso::TwoPartDate ut1 = mjd_utc;
  ut1._small += eops.dut / 86400e0; // add UT1-UTC, interpolated
  // const auto sofa_ut1jd = ut1.normalized().jd_sofa();
  era = iers2010::sofa::era00(/*sofa_ut1jd._big, sofa_ut1jd._small*/ut1.normalized());

  // Estimate s' [radians]
  const double sp = iers2010::sofa::sp00(/*sofajd._big, sofajd._small*/mjd_tt);

  // Form the polar motion matrix (W); note that we need angular units in
  // radians
  rpom =
      iers2010::sofa::pom00(dso::sec2rad(eops.xp), dso::sec2rad(eops.yp), sp);

  return 0;
}

// IAU 2006/2000A, CIO based, using X,Y series
int dso::gcrs2itrs(const dso::TwoPartDate &mjd_tai,
                   const dso::EopRecord &eops, double X, double Y, double s, double sprime,
                   Eigen::Matrix<double, 3, 3> &rc2i, double &era,
                   Eigen::Matrix<double, 3, 3> &rpom, double &xlod) noexcept {

  // TAI to TT
  // const dso::TwoPartDate mjd_tt(mjd_tai._big, mjd_tai._small + TAI2TTFDAYS);
  //const dso::TwoPartDate mjd_tt = mjd_tai.tai2tt();

  // TAI to UTC
  const dso::TwoPartDate mjd_utc = mjd_tai.tai2utc();

  // assign interpolated LOD value, [sec/day]
  xlod = eops.lod;

  // split date as in SOFA Date&Time idiom
  // auto sofajd = mjd_tt.jd_sofa();

  // X,Y coordinates of celestial intermediate pole from series based
  // on IAU 2006 precession and IAU 2000A nutation.
  //double X, Y;
  //iers2010::sofa::xy06(sofajd._big, sofajd._small, X, Y);

  // The CIO locator s, positioning the Celestial Intermediate Origin on
  // the equator of the Celestial Intermediate Pole, given the CIP's X,Y
  // coordinates. Compatible with IAU 2006/2000A precession-nutation.
  //const double s = iers2010::sofa::s06(sofajd._big, sofajd._small, X, Y);

  // Add CIP corrections (arcsec to radians)
  //X += dso::sec2rad(eops.dx);
  //Y += dso::sec2rad(eops.dy);

  // Form the celestial to intermediate-frame-of-date matrix given the CIP
  // X,Y and the CIO locator s.
  // [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]  = RC2T * [CRS]
  rc2i = iers2010::sofa::c2ixys(X, Y, s);

  // call ERA00 to get the ERA rotation angle (need UT1 datetime)
  dso::TwoPartDate ut1 = mjd_utc;
  ut1._small += eops.dut / 86400e0; // add UT1-UTC, interpolated
  //const auto sofa_ut1jd = ut1.normalized().jd_sofa();
  era = iers2010::sofa::era00(ut1.normalized());

  // Estimate s' [radians]
  const double sp = sprime;// iers2010::sofa::sp00(sofajd._big, sofajd._small);

  // Form the polar motion matrix (W); note that we need angular units in
  // radians
  rpom =
      iers2010::sofa::pom00(dso::sec2rad(eops.xp), dso::sec2rad(eops.yp), sp);

  return 0;
}
