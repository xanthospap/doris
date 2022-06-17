#include "astrodynamics.hpp"

Eigen::Matrix<double, 3, 1>
dso::solar_radiation_acceleration(const Eigen::Matrix<double, 3, 1> &rsat,
                                  const Eigen::Matrix<double, 3, 1> &rsun,
                                  double Area, double mass, double C) noexcept {
  // Relative position vector of spacecraft w.r.t. Sun
  const Eigen::Matrix<double, 3, 1> d = rsat - rsun;
  const double D = d.norm();
  const double D3 = D * D * D;
  constexpr const double AU2 = iers2010::AU * iers2010::AU;
  // Solar radiation pressure at 1 AU
  constexpr const double P0 =
      4.56316e-6; // [N/m^2] (~1367 W/m^2) from Beutler 2005

  // acceleration, Montenbruck 3.75
  return C * (Area / mass) * P0 * AU2 * d / D3;
}

Eigen::Matrix<double, 3, 1>
dso::solar_radiation_acceleration(const Eigen::Matrix<double, 3, 1> &rsat,
                                  const Eigen::Matrix<double, 3, 1> &rsun,
                                  double Area, double mass, double C,
                                  Eigen::Matrix<double, 3, 3> &dadr,
                                  Eigen::Matrix<double, 3, 1> &dadC) noexcept {
  // Relative position vector of spacecraft w.r.t. Sun
  const Eigen::Matrix<double, 3, 1> d = rsat - rsun;
  const double D = d.norm();
  const double D3 = D * D * D;
  constexpr const double AU2 = iers2010::AU * iers2010::AU;
  // Solar radiation pressure at 1 AU
  constexpr const double P0 =
      4.56316e-6; // [N/m^2] (~1367 W/m^2) from Beutler 2005

  // partials w.r.t sat position, see Monentbruck, 7.3.3
  dadr = P0 * C * (Area / mass) * AU2 *
         (Eigen::Matrix<double, 3, 3>::Identity() / D3 -
          3e0 * (d * d.transpose()) / (D3 * D * D));

  // acceleration
  const Eigen::Matrix<double, 3, 1> acceleration =
      C * (Area / mass) * P0 * AU2 * d / D3;

  // partials w.r.t C, see Monentbruck, 7.3.3
  dadC = acceleration / C;

  return acceleration;
}
