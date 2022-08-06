#include "orbit_integration.hpp"
#include "iers2010/iers2010.hpp"
#include "iers2010/iau.hpp"

Eigen::Matrix<double, 3, 3>
dso::itrs2gcrs(double mjd_tai, const dso::EopLookUpTable &eop_table,
          Eigen::Matrix<double, 3, 3> &ditrs2gcrs) noexcept {

  // split mjd to integral and fractional part
  double mjd_days;
  const double taif = std::modf(mjd_tai, &mjd_days);
  // get leap seconds
  const int leap_sec =
      dso::dat(dso::modified_julian_day(static_cast<int>(mjd_tai)));
  // compute fractional part of UTC date
  const double utcf = taif - static_cast<double>(leap_sec) / 86400e0;

  // interpolate/correct EOP values using UTC
  double xp, yp, dut1;
  if (eop_table.interpolate(mjd_days + utcf, xp, yp, dut1)) {
    fprintf(stderr, "ERROR. Failed getting EOP values");
  }

  // compute UT1 date (fractional part) UT1 = UTC + DUT1
  const double ut1f = utcf + (dut1 / 86400e0) * 1e-3;
  // compute TT date (fractional part) TT = TAI + 32.184sec
  const double ttf = taif + (32184e-3 / 86400e0);
  // UT1 date as MJD
  const double mjd_ut1 = mjd_days + ut1f;
  // TT date as MJD
  const double mjd_tt = mjd_days + ttf;

  // Form the celestial-to-intermediate matrix for this TT.
  auto rc2i = iers2010::sofa::c2i06a(dso::mjd0_jd, mjd_tt);
  
  // Predict the Earth rotation angle for this UT1.
  const double era = iers2010::sofa::era00(dso::mjd0_jd, mjd_ut1);
  
  // Estimate s'.
  const double sp = iers2010::sofa::sp00(dso::mjd0_jd, mjd_tt);
  
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
