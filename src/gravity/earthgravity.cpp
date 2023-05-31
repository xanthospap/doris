#include "egravity.hpp"
#include <stdexcept>

dso::EarthGravity::EarthGravity(const char *icgem, int degree, int order,
                                const dso::TwoPartDate &t)
    : max_degree(degree), max_order(order), cs_coeffs(degree),
      W(degree + 3, degree + 3), V(degree + 3, degree + 3) {
  if (dso::parse_gravity_model(icgem, degree, order, t, cs_coeffs)) {
    throw std::runtime_error(
        "[ERROR] Failed to construct Earth Gravity field\n");
  }
}

dso::iStatus dso::EarthGravity::acceleration(
    const Eigen::Matrix<double, 3, 1> &r_itrf,
    Eigen::Matrix<double, 3, 1> &acc_itrf,
    Eigen::Matrix<double, 3, 3> &gradient) noexcept {
  return dso::gravity_acceleration(cs_coeffs, r_itrf, max_degree,
                                   cs_coeffs.Re(), cs_coeffs.GM(), acc_itrf,
                                   gradient, &V, &W);
}
