#include "egravity.hpp"

Eigen::Matrix<double, 3, 1>
dso::point_mass_accel(double GM, const Eigen::Matrix<double, 3, 1> &rsat,
                      const Eigen::Matrix<double, 3, 1> &robj) noexcept {

  // Relative position vector of satellite w.r.t. point mass
  const Eigen::Matrix<double, 3, 1> d = rsat - robj;

  // we'll be needing some norms ...
  const double dn = d.norm();
  const double dn2 = d.squaredNorm();
  const double sn = robj.norm();
  const double sn2 = robj.squaredNorm();

  // acceleration
  return (-GM) * (d / (dn * dn2) + robj / (sn * sn2));
}

Eigen::Matrix<double, 3, 1>
dso::point_mass_accel(double GM, const Eigen::Matrix<double, 3, 1> &rsat,
                      const Eigen::Matrix<double, 3, 1> &robj,
                      Eigen::Matrix<double, 3, 3> &partials) noexcept {

  // Relative position vector of satellite w.r.t. point mass
  const Eigen::Matrix<double, 3, 1> d = rsat - robj;

  // we'll be needing some norms ...
  const double dn = d.norm();
  const double dn2 = d.squaredNorm();
  const double sn = robj.norm();
  const double sn2 = robj.squaredNorm();

  // acceleration
  // const Eigen::Matrix<double, 3, 1> a = (-GM) * (d / (dn*dn2) + robj /
  // (sn*sn2));

  // partials w.r.t. position
  const double sn5 = sn2 * sn2 * sn;
  partials = (-GM) * (Eigen::Matrix<double, 3, 3>::Identity() / (dn * dn2) +
                      3e0 * (d * d.transpose()) / sn5);

  return (-GM) * (d / (dn * dn2) + robj / (sn * sn2));
}
