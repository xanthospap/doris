#include "eop.hpp"
#include "orbit_integration.hpp"
#include "datetime/utcdates.hpp"
#include "iers2010/iers2010.hpp"
#include "iers2010/iau.hpp"
#include "eigen3/Eigen/Geometry"
#include <iers2010/iersc.hpp>

// TODO should have an error status!!!!
#ifdef ABCD
Eigen::Matrix<double, 3, 3>
dso::itrs2gcrs(double mjd_tai, const dso::EopLookUpTable &eop_table,
          Eigen::Matrix<double, 3, 3> &ditrs2gcrs) noexcept {

  // Need UTC datetime
  int imjd = (int)mjd_tai;
  double sec = (mjd_tai - (int)mjd_tai) * 86400e0;
  dso::nanoseconds::underlying_type iSec =
      static_cast<dso::nanoseconds::underlying_type>(
          sec * dso::nanoseconds::template sec_factor<double>());
  dso::datetime<dso::nanoseconds> taidate{dso::modified_julian_day(imjd),
                                          dso::nanoseconds(iSec)};
  dso::modified_julian_day utc_mjd;
  double utc = dso::tai2utc(taidate, utc_mjd);
  utc += static_cast<double>(utc_mjd.as_underlying_type());

  // interpolate/correct EOP values using UTC
  double xp, yp, dut1;
  if (int error; (error=eop_table.interpolate(utc, xp, yp, dut1))) {
    fprintf(stderr, "ERROR. Failed getting EOP values (status: %d)\n", error);
  }

  dso::datetime<dso::nanoseconds> ttdate = taidate;

  // Form the celestial-to-intermediate matrix for this TT.
  auto rc2i = iers2010::sofa::c2i06a(dso::mjd0_jd, ttdate.as_mjd());
  //printf("RC2I\n");
  //for (int i=0; i<3; i++) {
  //  for (int j=0; j<3; j++) {
  //    printf(" %+.6f ", rc2i(i,j));
  //  }
  //  printf("\n");
  //}
  
  // Predict the Earth rotation angle for this UT1.
  const double ut1 = utc + (dut1 / 86400e0) * 1e-3;
  const double era = iers2010::sofa::era00(dso::mjd0_jd, ut1);
  
  // Estimate s'.
  const double sp = iers2010::sofa::sp00(dso::mjd0_jd, ttdate.as_mjd());
  
  // Form the polar motion matrix.
  auto rpom =
      iers2010::sofa::pom00(xp * iers2010::DMAS2R, yp * iers2010::DMAS2R, sp);
  
  // Combine to form the celestial-to-terrestrial matrix.
  auto mat = rpom * dso::Mat3x3::RotZ(era) * rc2i;
  //printf("RC22\n");
  //for (int i=0; i<3; i++) {
  //  for (int j=0; j<3; j++) {
  //    printf(" %+.6f ", tmp(i,j));
  //  }
  //  printf("\n");
  //}
  
  // note that the following will result in an Eigen matrix that is the
  // transpose of mat (Eigen uses Column-Major and Mat3x3 Row-Major)
  Eigen::Matrix<double, 3, 3> t2c(mat.data);

  // ERA derivative
  const dso::Mat3x3 S({0e0, iers2010::OmegaEarth, 0e0, -iers2010::OmegaEarth,
                         0e0, 0e0, 0e0, 0e0, 0e0});
  mat = rpom * (S * dso::Mat3x3::RotZ(era)) * rc2i;
  ditrs2gcrs = Eigen::Matrix<double, 3, 3>(mat.data);

  return t2c;
}
#else

double OmegaEarth(double xlod) noexcept {
  return iers2010::OmegaEarth * (1e0 - xlod / 86400e0);
}

Eigen::Matrix<double, 3, 1>
dso::rcel2ter(const Eigen::Matrix<double, 3, 1> r,
                 const Eigen::Matrix<double, 3, 3> &rc2i,
            const double era,
                 const Eigen::Matrix<double, 3, 3> &rpom) noexcept {
  const auto rc2ti = Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) * rc2i;
  return (rpom * rc2ti) * r;
}

Eigen::Matrix<double, 6, 1>
dso::ycel2ter(const Eigen::Matrix<double, 6, 1> y,
              const Eigen::Matrix<double, 3, 3> &rc2i,
              const double era,
#ifdef NEW_EOP
              double xlod,
#endif
              const Eigen::Matrix<double, 3, 3> &rpom) noexcept {
  Eigen::Matrix<double, 6, 1> x;
  x.block<3, 1>(0, 0) = dso::rcel2ter(y.block<3, 1>(0, 0), rc2i, era, rpom);
  x.block<3, 1>(3, 0) = dso::rcel2ter(y.block<3, 1>(3, 0), rc2i, era, rpom);

  const double data[] = {0e0, -1e0, 0e0, 1e0, 0e0, 0e0, 0e0, 0e0, 0e0};

#ifdef NEW_EOP
  const auto rc2ti = (OmegaEarth(xlod) * Eigen::Matrix<double, 3, 3>(data) *
                      Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ())) *
                     rc2i;
#else
  const auto rc2ti = (iers2010::OmegaEarth * Eigen::Matrix<double, 3, 3>(data) *
                      Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ())) *
                     rc2i;
#endif
  x.block<3, 1>(3, 0) += (rpom * rc2ti) * y.block<3, 1>(0, 0);
  return x;
}

Eigen::Matrix<double, 3, 1>
dso::rter2cel(const Eigen::Matrix<double, 3, 1> r,
                 const Eigen::Matrix<double, 3, 3> &rc2i,
                 const double era,
                 const Eigen::Matrix<double, 3, 3> &rpom) noexcept {
  const auto rc2ti = Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) * rc2i;
  return (rpom * rc2ti).transpose() * r;
}

Eigen::Matrix<double, 6, 1>
dso::yter2cel(const Eigen::Matrix<double, 6, 1> y,
              const Eigen::Matrix<double, 3, 3> &rc2i,
              const double era,
#ifdef NEW_EOP
              double xlod,
#endif
              const Eigen::Matrix<double, 3, 3> &rpom) noexcept {
  Eigen::Matrix<double, 6, 1> x;
  x.block<3, 1>(0, 0) = dso::rter2cel(y.block<3, 1>(0, 0), rc2i, era, rpom);
  x.block<3, 1>(3, 0) = dso::rter2cel(y.block<3, 1>(3, 0), rc2i, era, rpom);

  const double data[] = {0e0, -1e0, 0e0, 1e0, 0e0, 0e0, 0e0, 0e0, 0e0};

#ifdef NEW_EOP
  const auto rc2ti = (OmegaEarth(xlod) * Eigen::Matrix<double, 3, 3>(data) *
                      Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ())) *
                     rc2i;
#else
  const auto rc2ti = (iers2010::OmegaEarth * Eigen::Matrix<double, 3, 3>(data) *
                      Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ())) *
                     rc2i;
#endif
  x.block<3, 1>(3, 0) += (rpom * rc2ti).transpose() * y.block<3, 1>(0, 0);
  return x;
}

// IAU 2006/2000A, CIO based, using X,Y series
int dso::gcrs2itrs(double mjd_tai, const dso::EopLookUpTable &eop_table,
                   Eigen::Matrix<double, 3, 3> &rc2i,
                   double &era,
                   Eigen::Matrix<double, 3, 3> &rpom, double &xlod) noexcept {

  // TAI MJD to datetime instance
  int imjd = (int)mjd_tai;
  double sec = (mjd_tai - (int)mjd_tai) * 86400e0;
  dso::nanoseconds::underlying_type iSec =
      static_cast<dso::nanoseconds::underlying_type>(
          sec * dso::nanoseconds::template sec_factor<double>());
  dso::datetime<dso::nanoseconds> taidate{dso::modified_julian_day(imjd),
                                          dso::nanoseconds(iSec)};

  // TAI to TT
  dso::datetime<dso::nanoseconds> ttdate(taidate);
  ttdate.add_seconds(dso::nanoseconds(32184 * 1'000'000L));

  // interpolate/correct EOP values using TT
  dso::EopRecord eops;
  if (int error; (error = eop_table.interpolate(ttdate.as_mjd(), eops))) {
    fprintf(stderr, "ERROR. Failed getting EOP values (status: %d)\n", error);
    return error;
  }

  // assign interpolated LOD value
  xlod = eops.lod;

  // X,Y coordinates of celestial intermediate pole from series based
  // on IAU 2006 precession and IAU 2000A nutation.
  double X, Y;
  iers2010::sofa::xy06(dso::mjd0_jd, ttdate.as_mjd(), X, Y);

  // The CIO locator s, positioning the Celestial Intermediate Origin on
  // the equator of the Celestial Intermediate Pole, given the CIP's X,Y
  // coordinates. Compatible with IAU 2006/2000A precession-nutation.
  const double s = iers2010::sofa::s06(dso::mjd0_jd, ttdate.as_mjd(), X, Y);

  // Add CIP corrections (arcsec to radians)
  X += dso::arcsec2rad(eops.dx);
  Y += dso::arcsec2rad(eops.dy);

  // Form the celestial to intermediate-frame-of-date matrix given the CIP
  // X,Y and the CIO locator s.
  // [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]  = RC2T * [CRS]
  rc2i = iers2010::sofa::c2ixys_e(X, Y, s);

  // call ERA00 to get the ERA rotation angle (need UT1 datetime)
  dso::modified_julian_day utc_mjd;
  const double utc = dso::tai2utc(taidate, utc_mjd); // fractional part
  double ut1 = utc + (eops.dut / 86400e0); // add UT1-UTC, interpolated
  ut1 += static_cast<double>(utc_mjd.as_underlying_type()); // UTC as mjd
  era = iers2010::sofa::era00(dso::mjd0_jd, ut1);

  // Estimate s' [radians]
  const double sp = iers2010::sofa::sp00(dso::mjd0_jd, ttdate.as_mjd());

  // Form the polar motion matrix (W); note that we need angular units in 
  // radians
  rpom = iers2010::sofa::pom00_e(dso::arcsec2rad(eops.xp),
                                 dso::arcsec2rad(eops.yp), sp);

  return 0;
}
#endif
