#include "orbit_integration.hpp"
#include "planetpos.hpp"

void dso::SunMoon(const dso::TwoPartDate &mjd_tai,
                  const Eigen::Matrix<double, 3, 1> &rsat, double GMSun,
                  double GMMoon, Eigen::Matrix<double, 3, 1> &sun_acc,
                  Eigen::Matrix<double, 3, 1> &moon_acc,
                  Eigen::Matrix<double, 3, 1> &sun_pos,
                  Eigen::Matrix<double, 3, 1> &mon_pos,
                  Eigen::Matrix<double, 3, 3> &gradient) noexcept {

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
  sun_acc = dso::point_mass_accel(GMSun, rsat, rSun, sgrad);

  /* Moon-induced acceleration [m/sec^2] and gradient */
  Eigen::Matrix<double, 3, 3> mgrad = Eigen::Matrix<double, 3, 3>::Zero();
  moon_acc = dso::point_mass_accel(GMMoon, rsat, rMon, mgrad);

  gradient = sgrad + mgrad;
  sun_pos = rSun;
  mon_pos = rMon;

  return;
}
