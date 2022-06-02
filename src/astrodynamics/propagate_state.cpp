#include "astrodynamics.hpp"
#include <matvec/matvec.hpp>
#include <cstdio>

using dso::Vector3;

int dso::propagate_state(double GM, const Vector3 &r0, const Vector3 &v0,
                         double dt, Vector3 &r, Vector3 &v) noexcept {

  // state to orbital elements
  dso::OrbitalElements ele;
  int error = state2elements(r0, v0, ele, GM);

  // propagate state to next epoch (dt)
  error = dso::elements2state(ele, dt, r, v, GM);

  return error;
}

Eigen::Matrix<double, 6, 1>
dso::propagate_state(double GM, const Eigen::Matrix<double, 6, 1> &Y0,
                     double dt) noexcept {
  const Vector3 r0({Y0(0), Y0(1), Y0(2)});
  const Vector3 v0({Y0(3), Y0(4), Y0(5)});
  Vector3 r, v;
  dso::propagate_state(GM, r0, v0, dt, r, v);
  Eigen::Matrix<double, 6, 1> Y;
  Y << r(0), r(1), r(2), v(0), v(1), v(2);
  return Y;
}

int dso::propagate_state(double GM, const Vector3 &r0, const Vector3 &v0,
                         double dt, Vector3 &r, Vector3 &v,
                         Eigen::Matrix<double, 6, 6> &dYdY0) noexcept {
  // state to orbital elements
  dso::OrbitalElements ele0;
  if (state2elements(r0, v0, ele0, GM))
    return 1;

  const double a = ele0.semimajor();
  const double e = ele0.eccentricity();
  const double i = ele0.inclination();

  // propagate state
  if (propagate_state(GM, r0, v0, dt, r, v))
    return 2;

  // State vector partials w.r.t epoch elements
  Eigen::Matrix<double, 6, 6> dY0dA0;
  if (dso::state_partials(ele0, GM, 0e0, dY0dA0))
    return 3;

  Eigen::Matrix<double, 6, 6> dYdA0;
  if (dso::state_partials(ele0, GM, dt, dYdA0))
    return 3;

  // cartesian to keplerian partial derivatives
  const double sqe2 = std::sqrt((1e0 - e) * (1e0 + e));
  const double n = std::sqrt(GM / (a * a * a));
  const double naa = std::sqrt(GM / a) * a;
  const double P_aM = -2e0 / (n * a); // P(a,M)     = -P(M,a)
  const double P_eM =
      -(1e0 - e) * (1e0 + e) / (naa * e); // P(e,M)     = -P(M,e)
  const double P_eo = +sqe2 / (naa * e);  // P(e,omega) = -P(omega,e)
  const double P_io =
      -1e0 / (naa * sqe2 * std::tan(i)); // P(i,omega) = -P(omega,i)
  const double P_iO =
      +1e0 / (naa * sqe2 * std::sin(i)); // P(i,Omega) = -P(Omega,i)

  // Partials of epoch elements w.r.t. epoch state
  Eigen::Matrix<double, 6, 6> dA0dY0;
  for (int j = 0; j < 3; j++) {
    dA0dY0(0, j) = P_aM * dY0dA0(j + 3, 5);
    dA0dY0(0, j + 3) = -P_aM * dY0dA0(j, 5);

    dA0dY0(1, j) = +P_eo * dY0dA0(j + 3, 4) + P_eM * dY0dA0(j + 3, 5);
    dA0dY0(1, j + 3) = -P_eo * dY0dA0(j, 4) - P_eM * dY0dA0(j, 5);

    dA0dY0(2, j) = +P_iO * dY0dA0(j + 3, 3) + P_io * dY0dA0(j + 3, 4);
    dA0dY0(2, j + 3) = -P_iO * dY0dA0(j, 3) - P_io * dY0dA0(j, 4);

    dA0dY0(3, j) = -P_iO * dY0dA0(j + 3, 2);
    dA0dY0(3, j + 3) = +P_iO * dY0dA0(j, 2);

    dA0dY0(4, j) = -P_eo * dY0dA0(j + 3, 1) - P_io * dY0dA0(j + 3, 2);
    dA0dY0(4, j + 3) = +P_eo * dY0dA0(j, 1) + P_io * dY0dA0(j, 2);

    dA0dY0(5, j) = -P_aM * dY0dA0(j + 3, 0) - P_eM * dY0dA0(j + 3, 1);
    dA0dY0(5, j + 3) = +P_aM * dY0dA0(j, 0) + P_eM * dY0dA0(j, 1);
  }

  // state transition matrix
  dYdY0 = dYdA0 * dA0dY0;

  return 0;
}

Eigen::Matrix<double, 6, 1>
dso::propagate_state(double GM, const Eigen::Matrix<double, 6, 1> &Y0,
                     double dt, Eigen::Matrix<double, 6, 6> &dYdY0) noexcept {
  const Vector3 r0({Y0(0), Y0(1), Y0(2)});
  const Vector3 v0({Y0(3), Y0(4), Y0(5)});
  Vector3 r, v;
  dso::propagate_state(GM, r0, v0, dt, r, v, dYdY0);
  Eigen::Matrix<double, 6, 1> Y;
  Y << r(0), r(1), r(2), v(0), v(1), v(2);
  return Y;
}
