#include "orbit_integration.hpp"
#include "datetime/utcdates.hpp"
#include "iers2010/iers2010.hpp"
#include "iers2010/iau.hpp"
#include <datetime/dtfund.hpp>

// TODO should have an error status!!!!

Eigen::Matrix<double, 3, 3>
dso::itrs2gcrs(double mjd_tai, const dso::EopLookUpTable &eop_table,
          Eigen::Matrix<double, 3, 3> &ditrs2gcrs) noexcept {

  // split mjd to integral and fractional part
  //double mjd_days;
  //const double taif = std::modf(mjd_tai, &mjd_days);
  //// get leap seconds
  //const int leap_sec =
  //    dso::dat(dso::modified_julian_day(static_cast<int>(mjd_tai)));
  //// compute fractional part of UTC date
  //const double utcf = taif - static_cast<double>(leap_sec) / 86400e0;

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
  if (int error; (error=eop_table.interpolate(/*mjd_days + utcf*/utc, xp, yp, dut1))) {
    fprintf(stderr, "ERROR. Failed getting EOP values (status: %d)\n", error);
  }

  //// compute UT1 date (fractional part) UT1 = UTC + DUT1
  //const double ut1f = utcf + (dut1 / 86400e0) * 1e-3;
  //// compute TT date (fractional part) TT = TAI + 32.184sec
  //const double ttf = taif + (32184e-3 / 86400e0);
  //// UT1 date as MJD
  //const double mjd_ut1 = mjd_days + ut1f;
  //// TT date as MJD
  //const double mjd_tt = mjd_days + ttf;
  dso::datetime<dso::nanoseconds> ttdate = taidate;

  // Form the celestial-to-intermediate matrix for this TT.
  auto rc2i = iers2010::sofa::c2i06a(dso::mjd0_jd, /*mjd_tt*/ttdate.as_mjd());
  
  // Predict the Earth rotation angle for this UT1.
  const double ut1 = utc + (dut1 / 86400e0) * 1e-3;
  const double era = iers2010::sofa::era00(dso::mjd0_jd, /*mjd_ut1*/ut1);
  
  // Estimate s'.
  const double sp = iers2010::sofa::sp00(dso::mjd0_jd, /*mjd_tt*/ttdate.as_mjd());
  
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
