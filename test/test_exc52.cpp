#include "astrodynamics.hpp"
#include "datetime/datetime_write.hpp"
#include "egravity.hpp"
#include "eigen3/Eigen/Eigen"
#include "geodesy/units.hpp"
#include "icgemio.hpp"
#include "iers2010/iau.hpp"
#include "iers2010/iers2010.hpp"
#include "iers_bulletin.hpp"
#include "integrators.hpp"
#include "planetpos.hpp"
#include "sp3/sp3.hpp"
#include <cassert>
#include <cstdio>
#include <datetime/dtfund.hpp>

const int dat = 13; //leap seconds at 1999/03/04[sec]a

using Datetime = dso::datetime<dso::nanoseconds>;

Datetime datetimeFromMjd(double mjd1, double mjd2 = 0e0) {
  double days;
  double fraction = std::modf(mjd1 + mjd2, &days);
  const double fdays2sec = fraction * 86400e0;
  const unsigned long nsec = static_cast<unsigned long>(fdays2sec * 1e9);
  return Datetime(dso::modified_julian_day((int)days), dso::nanoseconds(nsec));
}

dso::Mat3x3 RzMat(double angle) noexcept {
  const double s = std::sin(angle);
  const double c = std::cos(angle);
  return dso::Mat3x3({c, s, 0e0, -s, c, 0e0, 0e0, 0e0, 1e0});
}

Eigen::Matrix<double, 3, 3>
ter2cel(double mjd_gsp, Eigen::Matrix<double, 3, 3> *dt2c) noexcept {

  // keep small part, do computations with this
  double mjd_days;
  const double gpsf = std::modf(mjd_gsp, &mjd_days);
  const double taif = gpsf + (19e0 / 86400);
  const int leap_sec =
      dso::dat(dso::modified_julian_day(static_cast<int>(mjd_gsp)));
  if (leap_sec != dat) {
    fprintf(stderr, "Warning, found leap seconds to be: %d\n", leap_sec);
  }
  const double utcf = taif - static_cast<double>(leap_sec) / 86400e0;

  // we now have date in UTC, get EOP values
  double xp=0e0, yp=0e0, dut1=0e0;

  const double ut1f = utcf + (dut1 / 86400e0) * 1e-3;
  const double ttf = taif + (32184e-3 / 86400e0);

  // const double mjd_utc = mjd_days + utcf;
  const double mjd_ut1 = mjd_days + ut1f;
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
  auto mat = iers2010::sofa::c2tcio(rc2i, era, rpom);
  // note that the following will result in an Eigen matrix that is the
  // transpose of mat (Eigen uses Column-Major and Mat3x3 Row-Major)
  Eigen::Matrix<double, 3, 3> t2c(mat.data);

  /* ERA derivative */
  if (dt2c) {
    dso::Mat3x3 S({0e0, iers2010::OmegaEarth, 0e0, -iers2010::OmegaEarth,
                         0e0, 0e0, 0e0, 0e0, 0e0});
    mat = rpom * rc2i * S.rotz(era);
    *dt2c = Eigen::Matrix<double, 3, 3>(mat.data);
  }

  return t2c;
}

int main() {
  Datetime t1(dso::year(1999), dso::month(3), dso::day_of_month(4),
              dso::nanoseconds(0));

  Eigen::Matrix<double,3,3> dt2c;
  auto t2c = ter2cel(t1.as_mjd(), &dt2c);


  printf("U Matrix:\n");
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      printf(" %+.9f ", t2c(i,j));
    }
    printf("\n");
  }
  
  printf("dU Matrix: [10^4/sec]\n");
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      printf(" %+.9f ", dt2c(i,j)*1e4);
    }
    printf("\n");
  }
  
  return 0;
}
