#include "orbit_integration.hpp"
#include "planetpos.hpp"

void dso::SunMoon(const dso::TwoPartDate &mjd_tai,
                  const Eigen::Matrix<double, 3, 1> &rsat, double GMSun,
                  double GMMoon, Eigen::Matrix<double, 3, 1> &sun_acc,
                  Eigen::Matrix<double, 3, 1> &moon_acc,
                  Eigen::Matrix<double, 3, 1> &sun_pos,
                  Eigen::Matrix<double, 3, 1> &mon_pos,
                  Eigen::Matrix<double, 3, 3> &gradient) noexcept {

  /*
  const auto mjd_tt = mjd_tai.tai2tt();
  const double jd = (mjd_tt._big + dso::mjd0_jd) + mjd_tt._small;
  double rsun[3], rmon[3];

  // position vector of sun/moon, in J2000, [km]
  dso::cspice::j2planet_pos_from(dso::cspice::jd2et(jd), 10, 399, rsun);
  dso::cspice::j2planet_pos_from(dso::cspice::jd2et(jd), 301, 399, rmon);

  Eigen::Matrix<double, 3, 1> rSun(rsun); // [km]
  Eigen::Matrix<double, 3, 1> rMon(rmon); // [km]

  // Sun-induced acceleration [m/sec^2] and gradient
  Eigen::Matrix<double, 3, 3> sgrad = Eigen::Matrix<double, 3, 3>::Zero();
  sun_acc = dso::point_mass_accel(GMSun * 1e9, rsat, rSun * 1e3, sgrad);

  // Moon-induced acceleration [m/sec^2] and gradient
  Eigen::Matrix<double, 3, 3> mgrad = Eigen::Matrix<double, 3, 3>::Zero();
  moon_acc = dso::point_mass_accel(GMMoon * 1e9, rsat, rMon * 1e3, mgrad);

  // gradient
  gradient = sgrad + mgrad;

  // Sun position in [m]
  sun_pos = rSun * 1e3;
  mon_pos = rMon * 1e3;
  */

  /* Sun position in GCRF [m] */
  Eigen::Matrix<double, 3, 1> rSun;
  if (dso::planet_pos(dso::Planet::SUN, mjd_tai.tai2tt(), rSun)) {
    fprintf(stderr, "[ERROR] Failed getting sun position! (traceback: %s)\n",
            __func__);
    return;
  }
  /* Moon position in GCRF [m] */
  Eigen::Matrix<double, 3, 1> rMon; // [m]
  if (dso::planet_pos(dso::Planet::MOON, mjd_tai.tai2tt(), rMon)) {
    fprintf(stderr, "[ERROR] Failed getting moon position! (traceback: %s)\n",
            __func__);
    return;
  }

  /* Sun-induced acceleration [m/sec^2] and gradient */
  Eigen::Matrix<double, 3, 3> sgrad = Eigen::Matrix<double, 3, 3>::Zero();
  sun_acc = dso::point_mass_accel(GMSun * 1e9, rsat, rSun, sgrad);

  /* Moon-induced acceleration [m/sec^2] and gradient */
  Eigen::Matrix<double, 3, 3> mgrad = Eigen::Matrix<double, 3, 3>::Zero();
  moon_acc = dso::point_mass_accel(GMMoon * 1e9, rsat, rMon, mgrad);

  gradient = sgrad + mgrad;
  sun_pos = rSun;
  mon_pos = rMon;

  return;
}
