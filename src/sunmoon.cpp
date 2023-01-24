#include "orbit_integration.hpp"
#include "planetpos.hpp"
#include <cmath>
#include <datetime/dtcalendar.hpp>

void dso::SunMoon(/*double mjd_tai*/const dso::TwoPartDate &mjd_tai, const Eigen::Matrix<double, 3, 1> &rsat,
                  double GMSun, double GMMoon,
                  Eigen::Matrix<double, 3, 1> &sun_acc,
                  Eigen::Matrix<double, 3, 1> &moon_acc,
                  Eigen::Matrix<double, 3, 1> &sun_pos,
                  Eigen::Matrix<double, 3, 1> &mon_pos,
                  Eigen::Matrix<double, 3, 3> &mon_partials) noexcept {

  /*
  // split TAI to integral and fractional part
  double mjd_days;
  const double taif = std::modf(mjd_tai, &mjd_days);
  // compute (fractional part) of TT
  const double ttf = taif + (32184e-3 / 86400e0);
  // make TT date as MJD
  const double mjd_tt = mjd_days + ttf;
  // TT to MJD
  const double jd = mjd_tt + dso::mjd0_jd; // date as JD (TT)
  */
  const auto mjd_tt = mjd_tai.tai2tt();
  const double jd = (mjd_tt._big + dso::mjd0_jd) + mjd_tt._small;
  double rsun[3], rmon[3];

  // position vector of sun/moon, in J2000, [km]
  dso::cspice::j2planet_pos_from(dso::cspice::jd2et(jd), 10, 399, rsun);
  dso::cspice::j2planet_pos_from(dso::cspice::jd2et(jd), 301, 399, rmon);

  Eigen::Matrix<double, 3, 1> rSun(rsun); // [km]
  Eigen::Matrix<double, 3, 1> rMon(rmon); // [km]

  // Sun-induced acceleration [km/sec^2]
  sun_acc = dso::point_mass_accel(GMSun, rsat * 1e-3, rSun);
  sun_acc = sun_acc * 1e-3; // [m/sec^2]

  // Moon-induced acceleration [m/sec^2]
  moon_acc =
      dso::point_mass_accel(GMMoon * 1e9, rsat, rMon * 1e3, mon_partials);

  // Sun position in [m]
  sun_pos = rSun * 1e3;
  mon_pos = rMon * 1e3;

  return;
}
