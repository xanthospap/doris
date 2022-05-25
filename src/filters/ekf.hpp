#ifndef __EXTENDED_KALMAN_FILTER_HPP__
#define __EXTENDED_KALMAN_FILTER_HPP__

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Eigen"

namespace dso {

template <int N, typename S> struct ExtendedKalmanFilter {
  dso::datetime<S> t;
  Eigen::Matrix<double, N, 1> x;
  Eigen::Matrix<double, N, N> P;
  Eigen::Matrix<double, N, 1> K;

  ExtendedKalmanFilter(const dso::datetime<S> &t0,
                  const Eigen::Matrix<double, N, 1> &x0,
                  const Eigen::Matrix<double, N, N> &P0) noexcept {
    t = t0;
    x = x0;
    P = P0;
  }

  ExtendedKalmanFilter(const dso::datetime<S> &t0,
                  const Eigen::Matrix<double, N, 1> &x0,
                  const double *sigmas) noexcept {
    t = t0;
    x = x0;
    P.setZero();
    for (int i = 0; i < N; i++)
      P(i, i) = sigmas[i] * sigmas[i];
  }

  void time_update(const dso::datetime<S> &tk,
                   const Eigen::Matrix<double, N, 1> &xk,
                   const Eigen::Matrix<double, N, N> &phi) noexcept {
    t = tk;
    x = xk;
    P = phi * P * phi.transpose();
  }

  void time_update(const dso::datetime<S> &tk,
                   const Eigen::Matrix<double, N, 1> &xk,
                   const Eigen::Matrix<double, N, N> &phi,
                   const Eigen::Matrix<double, N, N> &Qdt) noexcept {
    t = tk;
    x = xk;
    P = phi * P * phi.transpose() + Qdt;
  }

  void observation_update(double z, double g, double sigma,
                          const Eigen::Matrix<double, N, 1> &H) noexcept {
    double inv_w = sigma * sigma;
    // kalman gain
    K = P * H / (inv_w + H * P * H);
    // state update
    x = x + K * (z - g);
    // covariance update (Joseph)
    auto KWm1Kt = (K * sigma) * (K * sigma).transpose();
    auto ImKG = Eigen::Matrix<double, N, N>::Identity() - K * H.transpose();
    P = ImKG * P * ImKG.transpose() + KWm1Kt;
  }

}; // ExtendedKalmanFilter

} // namespace dso

#endif
