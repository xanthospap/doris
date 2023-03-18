#ifndef __EXTENDED_KALMAN_FILTER_HPP__
#define __EXTENDED_KALMAN_FILTER_HPP__

#include "eigen3/Eigen/Eigen"
#include "datetime/dtcalendar.hpp"
#ifdef DEBUG
#include <cstdio>
#endif

namespace dso {

struct ExtendedKalmanFilter {
  dso::TwoPartDate t;
  Eigen::VectorXd x;
  Eigen::MatrixXd P;
  Eigen::VectorXd K;

  ExtendedKalmanFilter() noexcept {};

  ExtendedKalmanFilter(int num_params) noexcept;

  ExtendedKalmanFilter(const dso::TwoPartDate &t0, const Eigen::VectorXd &x0,
                       const Eigen::MatrixXd &P0) noexcept;

  const Eigen::VectorXd &state() const noexcept { return x; }

  dso::TwoPartDate time() const noexcept { return t; }

  void time_update(const dso::TwoPartDate &tk, const Eigen::VectorXd &xk,
                   const Eigen::MatrixXd &phi) noexcept;

  void time_update(const dso::TwoPartDate &tk, const Eigen::VectorXd &xk,
                   const Eigen::MatrixXd &phi,
                   const Eigen::MatrixXd &Qdt) noexcept;

  void observation_update(double z, double g, double sigma,
                          const Eigen::VectorXd &H) noexcept;

}; // ExtendedKalmanFilter

} // namespace dso

#endif
