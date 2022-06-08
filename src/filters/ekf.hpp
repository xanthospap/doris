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

  const Eigen::Matrix<double, N, 1> &state() const noexcept {return x;}
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
    //#ifdef DEBUG
    //printf("\tekf::time_update Matrix phi:\n");
    //for (int i=0;i<6;i++) {
    //  printf("\t\t|");
    //  for (int j=0; j<6; j++)
    //    printf("%8.1f ", phi(i,j));
    //  printf("\n");
    //}
    //#endif
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
    printf("\tObservation Update: z=%.3f g=%.3f sigma=%.3f\n", z,g,sigma);
    printf("\t                  : H=[%.3e %.3e %.3e %.3e %.3e %.3e]\n", H(0),H(1),H(2),H(3),H(4),H(5));
    printf("\t                  : P=\n");
    for (int k=0;k<6;k++) {
      printf("\t\t|");
      for (int j=0; j<6; j++)
        printf("%15.5f ", P(k,j));
      printf("\n");
    }
    double inv_w = sigma * sigma;
    // kalman gain
    // K = P * H / (inv_w + H.transpose() * P * H);
    K = P * H / (inv_w + H.dot(P * H));
    const auto PH = P*H;
    printf("\t                  :PH=[%.3f %.3f %.3f %.3f %.3f %.3f]\n", PH(0),PH(1),PH(2),PH(3),PH(4),PH(5));
    printf("\t                  : K=[%.3f %.3f %.3f %.3f %.3f %.3f]\n", K(0),K(1),K(2),K(3),K(4),K(5));
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
