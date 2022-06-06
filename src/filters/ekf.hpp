#ifndef __EXTENDED_KALMAN_FILTER_HPP__
#define __EXTENDED_KALMAN_FILTER_HPP__

#include "iers2010/matvec.hpp"
#include "eigen3/Eigen/Eigen"
#ifdef DEBUG
#include <cstdio>
#endif

namespace dso {

template <int N, typename S> struct ExtendedKalmanFilter {
  dso::datetime<S> t;
  Eigen::Matrix<double, N, 1> x;
  Eigen::Matrix<double, N, N> P;
  Eigen::Matrix<double, N, 1> K;

  Vector3 state_position_vector() const noexcept {
    return Vector3({x(0), x(1), x(2)});
  }
  
  Vector3 state_velocity_vector() const noexcept {
    return Vector3({x(3), x(4), x(5)});
  }

  Eigen::Matrix<double, N, 1> state() const noexcept {return x;}
  dso::datetime<S> time() const noexcept {return t;}

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
    #ifdef DEBUG
      //printf("\tEKF::time update ->\n");
      //printf("\tx=[%+12.2f %+12.2f %+12.2f %+12.2f %+12.2f %+12.2f]\n", x(0),x(1),x(2),x(3),x(4),x(5));
      //printf("\t  [%+6.3f %+6.3f %+6.3f %+6.3f %+6.3f %+6.3f]\n", phi(0,0),phi(0,1),phi(0,2),phi(0,3),phi(0,4),phi(0,5));
      //printf("\t  [             %+6.3f %+6.3f %+6.3f %+6.3f %+6.3f]\n", phi(1,1),phi(1,2),phi(1,3),phi(1,4),phi(1,5));
      //printf("\t  [                          %+6.3f %+6.3f %+6.3f %+6.3f]\n", phi(2,2),phi(2,3),phi(2,4),phi(2,5));
      //printf("\t  [                                       %+6.3f %+6.3f %+6.3f]\n", phi(3,3),phi(3,4),phi(3,5));
      //printf("\t  [                                                    %+6.3f %+6.3f]\n", phi(4,4),phi(4,5));
      //printf("\t  [                                                                 %+6.3f]\n", phi(5,5));
    #endif
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
    // K = P * H / (inv_w + H.transpose() * P * H);
    K = P * H / (inv_w + H.dot(P * H));
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
