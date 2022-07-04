#include "rkf45.hpp"
#include <cmath>
#include <limits>
#ifdef DEBUG
#include <cassert>
#endif

/// @param[in] t0  Initial time (aka t0)
/// @param[in] y0  Initial values at time t
/// @param[in] yp0 Initial derivatives at time t
/// @param[in] h    Step size
/// @returns The fifth order (sixth order accurate locally) solution
///          approximation at t+h
Eigen::VectorXd dso::RKF45::fehl(double t0, double hstep,
                                 const Eigen::VectorXd &y0,
                                 const Eigen::VectorXd &yp0) noexcept {
  Eigen::VectorXd s(neqn);
#ifdef DEBUG
  assert(y0.rows() == neqn);
  assert(y0.cols() == 1);
  assert(yp0.rows() == neqn);
  assert(yp0.cols() == 1);
  assert(F.col(0).rows() == neqn);
#endif

  double ch = hstep / 4e0;
  F.col(4) = y0 + ch * yp0;
  f(t0 + ch, F.col(4), F.col(0), params);

  ch = 3e0 * hstep / 32e0;
  F.col(4) = y0 + ch * (yp + 3e0 * F.col(0));
  f(t0 + 3e0 * hstep / 8e0, F.col(4), F.col(1), params);

  ch = hstep / 2197e0;
  F.col(4) = y0 + ch * (1932e0 * yp + (7296e0 * F.col(1) - 7200e0 * F.col(0)));
  f(t0 + 12e0 * hstep / 13e0, F.col(4), F.col(2), params);

  ch = hstep / 4104e0;
  F.col(4) = y0 + ch * ((8341e0 * yp - 845e0 * F.col(2)) +
                       (29440e0 * F.col(1) - 32832e0 * F.col(0)));
  f(t0 + hstep, F.col(4), F.col(3), params);

  ch = hstep / 20520e0;
  F.col(0) =
      y0 + ch * ((-6080e0 * yp + (9295e0 * F.col(2) - 5643e0 * F.col(3))) +
                (41040e0 * F.col(0) - 28352e0 * F.col(1)));
  f(t0 + hstep / 2e0, F.col(0), F.col(4), params);

  // compute approximate solution at t+h

  ch = hstep / 7618050e0;
  return y0 +
         ch * ((902880e0 * yp + (3855735e0 * F.col(2) - 1371249e0 * F.col(3))) +
               (3953664e0 * F.col(1) + 277020e0 * F.col(4)));
}
