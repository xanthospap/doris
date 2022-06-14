#include "egravity.hpp"

Eigen::Matrix<double, 3, 1>
dso::point_mass_accel(double GM, const Eigen::Matrix<double, 3, 1> &rsat,
                      const Eigen::Matrix<double, 3, 1> &robj,
                      Eigen::Matrix<double, 3, 3> &partials) noexcept {
  const Eigen::Matrix<double, 3, 1> d = rsat - robj;
  const double D2 = d.dot(d);
  const double D = std::sqrt(D2);
  const double D3 = D2 * D;
  const double S2 = robj.dot(robj);
  const double S = std::sqrt(S2);
  const double S3 = S2 * S;
  const Eigen::Matrix<double, 3, 1> a = -GM * (d / D3 + robj / S3);
  const double D5 = D3 * D2;
  partials = -GM * (Eigen::Matrix<double, 3, 3>::Identity() / D3 -
                    3e0 * (d*d.transpose()) / D5);
  return a;
}
