#include "astrodynamics.hpp"
#include "gcem.hpp" // for constexpr math (trigonometric funcs)
#include "geodesy/ellipsoid.hpp"
#include <cmath>
#ifdef DEBUG
#include <cstdio>
#endif

//[[maybe_unused]]
// constexpr const double MeanRe =
// dso::mean_earth_radius<dso::ellipsoid::grs80>();

/// @file occultation.cpp
/// @brief Satellite eclipse, shadow and occultation models
/// @todo At the moment of writing, trigonometric std functions are not
///       constexpr. When they are, drop dependendy of gcem

double utest::vallado_shadow(const Eigen::Matrix<double,3,1> &r_sat,
                             const Eigen::Matrix<double,3,1> &r_sun) noexcept {
  const double product = r_sat.dot(r_sun);
  if (product >= 0)
    return 1e0;

  const double r = r_sat.norm();
  const double s = r_sun.norm();
  const double cosz = product / (r * s);
  const double zeta = std::acos(cosz);
  const double sinz = std::sin(zeta);
  const double sath = r * cosz;
  const double satv = r * sinz;
  constexpr const double sinApen = (iers2010::Rs + iers2010::Re) / iers2010::AU;
  constexpr const double Apen = /*std::asin(sinApen);*/ gcem::asin(sinApen);
  constexpr const double tanApen = /*std::tan(Apen);*/ gcem::tan(Apen);
  constexpr const double x = iers2010::Re / sinApen;
  const double penv = tanApen * (x + sath);
  if (satv <= penv) {
    // penumbra ...
    constexpr const double sinAumb =
        (iers2010::Rs - iers2010::Re) / iers2010::AU;
    constexpr const double y = iers2010::Re / sinAumb;
    constexpr const double Aumb = /*std::asin(sinAumb);*/ gcem::asin(sinAumb);
    constexpr const double tanAumb = /*std::tan(Aumb);*/ gcem::tan(Aumb);
    const double umbv = tanAumb * (y - sath);
    return (satv <= umbv) ? 0 : 0.5e0;
  }
  return 1;
}

/*
double utest::montebruck_shadow(const Eigen::Matrix<double,3,1> &r_sat,
                                const Eigen::Matrix<double,3,1> &r_sun) noexcept {
  auto unitvec_sun = r_sun / r_sun.norm();
  double s = r_sat.dot(unitvec_sun);
  return ((s > 0e0 || (r_sat - s * unitvec_sun).norm() > iers2010::Re) ? 1e0
                                                                       : 0e0);
}
*/
double
utest::montebruck_shadow(const Eigen::Matrix<double, 3, 1> &r_sat,
                         const Eigen::Matrix<double, 3, 1> &r_sun) noexcept {
  const auto unitvec_sun = r_sun.normalized(); // r_sun / r_sun.norm();
  const double s = r_sat.dot(unitvec_sun); // r_sat.dot_product(unitvec_sun);
  return ((s > 0e0 || (r_sat - s * unitvec_sun).norm() > iers2010::Re) ? 1e0
                                                                       : 0e0);
}

double utest::bernese_shadow1(const Eigen::Matrix<double,3,1> &r_sat,
                              const Eigen::Matrix<double,3,1> &r_sun) noexcept {
  constexpr const double fac =
      1e0 + dso::ellipsoid_traits<dso::ellipsoid::grs80>::f;
  Eigen::Matrix<double,3,1> zsat;zsat<<r_sat.x(), r_sat.y(), r_sat.z() * fac;
  Eigen::Matrix<double,3,1> zsun;zsun<<r_sun.x(), r_sun.y(), r_sun.z() * fac;
  const double rsun = zsun.norm();
  const double rsat2 = zsat.squaredNorm();
  const double rcos = zsat.dot(zsun) / rsun;
  const double rsin = std::sqrt(rsat2 - rcos * rcos);
  if (rcos < 0e0 && (rsin - iers2010::Re < 0e0))
    return 0e0;
  return 1;
}

double utest::conical_shadow(const Eigen::Matrix<double,3,1> &r_sat,
                             const Eigen::Matrix<double,3,1> &r_sun) noexcept {
  const double te = std::asin(iers2010::Re / r_sat.norm());
  const double nsunsat = (r_sun - r_sat).norm();
  const double ts = std::asin(iers2010::Rs / nsunsat);
  const double tes =
      std::acos(-r_sat.dot(r_sun - r_sat) / (nsunsat * r_sat.norm()));

  if (tes >= te + ts)
    return 1e0;
  if (tes <= te - ts)
    return 0e0;

  const double ts2 = ts * ts;
  const double te2 = te * te;
  const double tes2 = tes * tes;
  const double caf = std::acos((ts2 + tes2 - te2) / (2e0 * ts * tes));
  const double cbd = std::acos((te2 + tes2 - ts2) / (2e0 * te * tes));
  const double Safc = 0.5 * caf * ts2;
  const double Saec = 0.5 * (ts * std::sin(caf)) * (ts * std::cos(caf));
  const double Sbdc = 0.5 * cbd * te2;
  const double Sbec = 0.5 * (te * std::sin(cbd)) * (te * std::cos(cbd));
  const double S = 2e0 * (Safc - Saec) + 2e0 * (Sbdc - Sbec);
  return 1e0 - S / (iers2010::DPI * ts2);
}
