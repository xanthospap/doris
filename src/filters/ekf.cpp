#include "ekf.hpp"

dso::ExtendedKalmanFilter::ExtendedKalmanFilter(int num_params) noexcept
    : t(dso::TwoPartDate{}), x(Eigen::VectorXd::Zero(num_params)),
      P(Eigen::MatrixXd::Identity(num_params, num_params)),
      K(Eigen::VectorXd::Zero(num_params)) {}

dso::ExtendedKalmanFilter::ExtendedKalmanFilter(
    const dso::TwoPartDate &t0, const Eigen::VectorXd &x0,
    const Eigen::MatrixXd &P0) noexcept
    : t(t0), x(x0), P(P0), K(Eigen::VectorXd(x0.cols())) {}

void dso::ExtendedKalmanFilter::time_update(
    const dso::TwoPartDate &tk, const Eigen::VectorXd &xk,
    const Eigen::MatrixXd &phi) noexcept {
  t = tk;
  x = xk;
  P = phi * P * phi.transpose();
}

void dso::ExtendedKalmanFilter::time_update(
    const dso::TwoPartDate &tk, const Eigen::VectorXd &xk,
    const Eigen::MatrixXd &phi, const Eigen::MatrixXd &Qdt) noexcept {
  t = tk;
  x = xk;
  P = phi * P * phi.transpose() + Qdt;
}

void dso::ExtendedKalmanFilter::observation_update(
    double z, double g, double sigma, const Eigen::VectorXd &H) noexcept {
  double inv_w = sigma * sigma;

  // kalman gain
  // K = P * H / (inv_w + H.transpose() * P * H);
  K = P * H / (inv_w + H.dot(P * H));

  // state update
  x = x + K * (z - g);

  const int n = x.rows();

  // covariance update (Joseph variant)
  auto KWm1Kt = (K * sigma) * (K * sigma).transpose();
  auto ImKG = Eigen::MatrixXd::Identity(n, n) - K * H.transpose();
  P = ImKG * P * ImKG.transpose() + KWm1Kt;
}
