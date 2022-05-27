#include "astrodynamics.hpp"
#include <iers2010/matvec.hpp>
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

int dso::propagate_state(double GM, const Vector3 &r0, const Vector3 &v0,
                         double dt, Vector3 &r, Vector3 &v,
                         Eigen::Matrix<double, 6, 6> &dYdY0) noexcept {
  // state to orbital elements
  dso::OrbitalElements ele0;
  if (state2elements(r0, v0, ele0, GM)) return 1;

  const double a = ele0.semimajor();
  const double e = ele0.eccentricity();
  const double i = ele0.inclination();

  // propagate state
  if (propagate_state(GM, r0, v0, dt, r, v)) return 2;

  // State vector partials w.r.t epoch elements
  Eigen::Matrix<double, 6, 6> dY0dA0;
  if (dso::state_partials(ele0, GM, 0e0, dY0dA0)) return 3;

  Eigen::Matrix<double, 6, 6> dYdA0;
  if (dso::state_partials(ele0, GM, dt, dYdA0)) return 3;

  // cartesian to keplerian partial derivatives
  const double sqe2 = std::sqrt((1e0-e)*(1e0+e));
  const double naa  = std::sqrt(GM/a)*a;
  const double P_aM = -2e0/(n*a);                   // P(a,M)     = -P(M,a)
  const double P_eM = -(1e0-e)*(1e0+e)/(naa*e);     // P(e,M)     = -P(M,e)
  const double P_eo = +sqe2/(naa*e);                // P(e,omega) = -P(omega,e)
  const double P_io = -1e0/(naa*sqe2*std::tan(i));  // P(i,omega) = -P(omega,i)
  const double P_iO = +1e0/(naa*sqe2*std::sin(i));  // P(i,Omega) = -P(Omega,i)


  // Partials of epoch elements w.r.t. epoch state
  for (int i=0; i<3; i++) {
    dA0dY0(0,i) = P_am * dY0dA0(i+3,5);
    dA0dY0(0,i+3) = - P_aM*dY0dA0(i, 5);

    dA0dY0(1,i)   = + P_eo*dY0dA0(i+3,4) + P_eM*dY0dA0(i+3,5);
    dA0dY0(1,i+3) = - P_eo*dY0dA0(i  ,4) - P_eM*dY0dA0(i  ,5);

    dA0dY0(2,i)   = + P_iO*dY0dA0(i+3,3) + P_io*dY0dA0(i+3,4);
    dA0dY0(2,i+3) = - P_iO*dY0dA0(i  ,3) - P_io*dY0dA0(i  ,4);

    dA0dY0(3,i)   = - P_iO*dY0dA0(i+3,2);
    dA0dY0(3,i+3) = + P_iO*dY0dA0(i  ,2);

    dA0dY0(4,i)   = - P_eo*dY0dA0(i+3,1) - P_io*dY0dA0(i+3,2);
    dA0dY0(4,i+3) = + P_eo*dY0dA0(i  ,1) + P_io*dY0dA0(i  ,2);

    dA0dY0(5,i)   = - P_aM*dY0dA0(i+3,0) - P_eM*dY0dA0(i+3,1);
    dA0dY0(5,i+3) = + P_aM*dY0dA0(i  ,0) + P_eM*dY0dA0(i  ,1);
  }

  // state transition matrix
  dYdY0 = dYdA0 * dA0dY0;
}
