#include "egravity.hpp"
#include "iers2010/iersc.hpp"
#include <stdexcept>

dso::EarthGravity::EarthGravity(const char *icgem, int degree, int order,
                                const dso::TwoPartDate &t)
    : cs_coeffs(degree, order, iers2010::GMe, iers2010::Re) {
  if (dso::parse_gravity_model(icgem, degree, order, t, cs_coeffs)) {
    throw std::runtime_error(
        "[ERROR] Failed to construct Earth Gravity field\n");
  }
}

dso::iStatus dso::EarthGravity::acceleration(
    const Eigen::Matrix<double, 3, 1> &r_itrf,
    Eigen::Matrix<double, 3, 1> &acc_itrf,
    Eigen::Matrix<double, 3, 3> &gradient, int max_degree, int max_order,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *V,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *W) noexcept {

  /* set max degree and order */
  if (max_degree < 0)
    max_degree = this->max_degree();
  if (max_order < 0)
    max_order = this->max_order();

  /* validate degree and order */
  if (max_degree > this->max_degree() || max_order > this->max_order() ||
      max_degree < max_order) {
    fprintf(stderr,
            "[ERROR] Invalid degree/order for gravity field SH computation "
            "(traceback: %s)\n",
            __func__);
    return dso::iStatus(1);
  }

  /* setup and validate workspace (if needed) */
  if (!V) {
    V = new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
        max_degree + 3, max_degree + 3);
  }
  if (!W) {
    W = new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
        max_degree + 3, max_degree + 3);
  }

  /* check sizes of workspace matrices */
  if (V->rows() < max_degree + 2 || W->rows() < max_degree + 2) {
    fprintf(stderr,
            "[ERROR] Invalid workspace matrices for SH computation "
            "(traceback: %s)\n",
            __func__);
    return dso::iStatus(1);
  }

  return dso::gravity_acceleration(cs_coeffs, r_itrf, max_degree,
                                   cs_coeffs.Re(), cs_coeffs.GM(), acc_itrf,
                                   gradient, V, W);
}
