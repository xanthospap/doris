#include "orbit_integration.hpp"
#include "datetime/utcdates.hpp"
#include "iers2010/iers2010.hpp"
#include "iers2010/iau.hpp"
#include "eigen3/Eigen/Geometry"

// TODO should have an error status!!!!

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

  
Eigen::Matrix<double, 6, 1>
dso::itrs2gcrs(const Eigen::Matrix<double, 6, 1> &y_itrs, double mjd_tai,
          const dso::EopLookUpTable &eop_table) noexcept {
  
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

  // call XY06 to get X, Y (series)
  double X, Y;
  iers2010::sofa::xy06(dso::mjd0_jd, ttdate.as_mjd(), X, Y);

  // call S06 to get s
  const double s = iers2010::sofa::s06(dso::mjd0_jd, ttdate.as_mjd(), X, Y);

  // TODO need to correct the X, Y values (see IERS1010, 5.5.4
  // Forced motion of the Celestial Intermediate Pole in the GCRS)

  // call C2IXYS to get the GCRS-to-CIRS matrix
  const Eigen::Matrix<double,3,3> gcrs2cirs = iers2010::sofa::c2ixys_e(X, Y, s);

  // Need UTC datetime
  dso::modified_julian_day utc_mjd;
  double utc = dso::tai2utc(taidate, utc_mjd);
  utc += static_cast<double>(utc_mjd.as_underlying_type());

  // interpolate/correct EOP values using UTC
  double xp, yp, dut1;
  if (int error; (error=eop_table.interpolate(utc, xp, yp, dut1))) {
    fprintf(stderr, "ERROR. Failed getting EOP values (status: %d)\n", error);
  }

  // call ERA00 to get the ERA rotation angle
  const double ut1 = utc + (dut1 / 86400e0) * 1e-3;
  const double era = iers2010::sofa::era00(dso::mjd0_jd, ut1);

  // dso::Mat3x3::RotZ(era) * rc2i;
  const Eigen::Matrix<double, 3, 3> pnr =
      Eigen::AngleAxisd(-era, Eigen::Vector3d::UnitZ()) * gcrs2cirs;

  // Estimate s'
  const double sp = iers2010::sofa::sp00(dso::mjd0_jd, ttdate.as_mjd());
  
  // Form the polar motion matrix (W)
  const Eigen::Matrix<double, 3, 3> rpom =
      iers2010::sofa::pom00_e(xp * iers2010::DMAS2R, yp * iers2010::DMAS2R, sp);

  // split state and position and velocity vectors
  const Eigen::Matrix<double,3,1> r =  y_itrs.block<3,1>(0,0);
  const Eigen::Matrix<double,3,1> v =  y_itrs.block<3,1>(3,0);
  const double _od[] = {0e0, 0e0, iers2010::OmegaEarth};
  const Eigen::Matrix<double,3,1> omega(_od);
  
  // 
  Eigen::Matrix<double,6,1> y_gcrs;
  y_gcrs.block<3,1>(0,0) = gcrs2cirs * rpom * r;
  y_gcrs.block<3,1>(3,0) = gcrs2cirs * (rpom * v + omega.cross(r));

  return y_gcrs;
}
