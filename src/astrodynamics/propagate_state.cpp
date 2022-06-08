#include "astrodynamics.hpp"
// #include <matvec/matvec.hpp>
#include <cstdio>

// using dso::Vector3;
/*
int dso::propagate_state(double GM, const Vector3 &r0, const Vector3 &v0,
                         double dt, Vector3 &r, Vector3 &v) noexcept {

  // state to orbital elements
  dso::OrbitalElements ele;
  int error = state2elements(r0, v0, ele, GM);

  // propagate state to next epoch (dt)
  error = dso::elements2state(ele, dt, r, v, GM);

  return error;
}
*/
Eigen::Matrix<double, 6, 1>
dso::propagate_state(double GM, const dso::OrbitalElements &Ele,
                     double dt) noexcept {
  // elements for ease of use
  const double a = Ele.semimajor();
  const double e = Ele.eccentricity();
  const double i = Ele.inclination();
  const double Omega = Ele.Omega();
  const double omega = Ele.omega();
  const double M0 = Ele.mean_anomaly();

  // compute mean anomaly M
  const double n = std::sqrt(GM / (a * a * a));
  const double M = M0 + n * dt;

  // eccentric anomaly
  int error;
  const double E = dso::kepler(e, M, error);
  const double cE = std::cos(E);
  const double sE = std::sin(E);

  // perifocal coordinates
  const double fac = std::sqrt((1e0 - e) * (1e0 + e));
  const double R = a * (1e0 - e * cE);    // Distance
  const double V = std::sqrt(GM * a) / R; // Velocity
  Eigen::Matrix<double, 6, 1> state;
  state << a * (cE - e), a * fac * sE, 0e0, -V * sE, +V * fac * cE, 0e0;

  // Transformation to reference system (Gaussian vectors)
  // (note we are using positive angle values here - no need to take the
  // transpose)
  Eigen::Matrix3d PQW(Eigen::AngleAxisd(Omega, Eigen::Vector3d::UnitZ()) *
                      Eigen::AngleAxisd(i, Eigen::Vector3d::UnitX()) *
                      Eigen::AngleAxisd(omega, Eigen::Vector3d::UnitZ()));

  state.block<3, 1>(0, 0) = PQW * state.block<3, 1>(0, 0);
  state.block<3, 1>(3, 0) = PQW * state.block<3, 1>(3, 0);

  return state;
}
/*
int dso::propagate_state(double GM, const Vector3 &r0, const Vector3 &v0,
                         double dt, Vector3 &r, Vector3 &v,
                         Eigen::Matrix<double, 6, 6> &dYdY0) noexcept {
  int status = 0;

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
  Eigen::Matrix<double, 6, 6> dY0dA0 = dso::state_partials(ele0, GM, 0e0, status);

  printf("\t%s -> dY0dA0\n\t", __func__);
  for (int k=0; k<6; k++) printf("%+6.1f ", dY0dA0(0,k));
  printf("\n");

  Eigen::Matrix<double, 6, 6> dYdA0 = dso::state_partials(ele0, GM, dt, status);

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
  for (int k = 0; k < 3; k++) {
        dA0dY0(0,k)   = + P_aM*dY0dA0(k+3,5);
    dA0dY0(0,k+3) = - P_aM*dY0dA0(k  ,5);

    dA0dY0(1,k)   = + P_eo*dY0dA0(k+3,4) + P_eM*dY0dA0(k+3,5);
    dA0dY0(1,k+3) = - P_eo*dY0dA0(k  ,4) - P_eM*dY0dA0(k  ,5);

    dA0dY0(2,k)   = + P_iO*dY0dA0(k+3,3) + P_io*dY0dA0(k+3,4);
    dA0dY0(2,k+3) = - P_iO*dY0dA0(k  ,3) - P_io*dY0dA0(k  ,4);

    dA0dY0(3,k)   = - P_iO*dY0dA0(k+3,2);
    dA0dY0(3,k+3) = + P_iO*dY0dA0(k  ,2);

    dA0dY0(4,k)   = - P_eo*dY0dA0(k+3,1) - P_io*dY0dA0(k+3,2);
    dA0dY0(4,k+3) = + P_eo*dY0dA0(k  ,1) + P_io*dY0dA0(k  ,2);

    dA0dY0(5,k)   = - P_aM*dY0dA0(k+3,0) - P_eM*dY0dA0(k+3,1);
    dA0dY0(5,k+3) = + P_aM*dY0dA0(k  ,0) + P_eM*dY0dA0(k  ,1);
  }

  // state transition matrix
  dYdY0 = dYdA0 * dA0dY0;

  return 0;
}
*/

Eigen::Matrix<double, 6, 1>
dso::propagate_state(double GM, const Eigen::Matrix<double, 6, 1> &Y0,
                     double dt, Eigen::Matrix<double, 6, 6> &dYdY0) noexcept {
  // state to orbital elements
  const auto elements = dso::state2elements(GM, Y0);

  // propagate state
  // const auto Y = dso::propagate_state(GM, elements, dt);

  // State vector partials w.r.t epoch elements
  const auto dY0dA0 = dso::state_partials(GM,elements,0e0);
  const auto dYdA0  = dso::state_partials(GM,elements,dt);
  //  printf("\tdso::propagate_state Matrix dYdA0: (note dt=%.3f)\n",dt);
  //  for (int i=0;i<6;i++) {
  //    printf("\t\t|");
  //    for (int j=0; j<6; j++)
  //      printf("%8.1f ", dYdA0(i,j));
  //    printf("\n");
  //  }

  const double a = elements.semimajor();
  const double e = elements.eccentricity();
  const double i = elements.inclination();
  const double n = std::sqrt(GM/(a*a*a));

  // Poisson brackets
  const double sqe2 = std::sqrt((1e0-e)*(1e0+e));
  const double naa  = n*a*a;

  const double P_aM = -2e0/(n*a);                   // P(a,M)     = -P(M,a)
  const double P_eM = -(1e0-e)*(1e0+e)/(naa*e);     // P(e,M)     = -P(M,e)
  const double P_eo = +sqe2/(naa*e);                // P(e,omega) = -P(omega,e)
  const double P_io = -1e0/(naa*sqe2*std::tan(i));  // P(i,omega) = -P(omega,i)
  const double P_iO = +1e0/(naa*sqe2*std::sin(i));  // P(i,Omega) = -P(Omega,i)

  // Partials of epoch elements w.r.t. epoch state
  Eigen::Matrix<double,6,6> dA0dY0;
  for (int k=0;k<3;k++) {
    dA0dY0(0,k)   = + P_aM*dY0dA0(k+3,5);
    dA0dY0(0,k+3) = - P_aM*dY0dA0(k  ,5);

    dA0dY0(1,k)   = + P_eo*dY0dA0(k+3,4) + P_eM*dY0dA0(k+3,5);
    dA0dY0(1,k+3) = - P_eo*dY0dA0(k  ,4) - P_eM*dY0dA0(k  ,5);

    dA0dY0(2,k)   = + P_iO*dY0dA0(k+3,3) + P_io*dY0dA0(k+3,4);
    dA0dY0(2,k+3) = - P_iO*dY0dA0(k  ,3) - P_io*dY0dA0(k  ,4);

    dA0dY0(3,k)   = - P_iO*dY0dA0(k+3,2);
    dA0dY0(3,k+3) = + P_iO*dY0dA0(k  ,2);

    dA0dY0(4,k)   = - P_eo*dY0dA0(k+3,1) - P_io*dY0dA0(k+3,2);
    dA0dY0(4,k+3) = + P_eo*dY0dA0(k  ,1) + P_io*dY0dA0(k  ,2);

    dA0dY0(5,k)   = - P_aM*dY0dA0(k+3,0) - P_eM*dY0dA0(k+3,1);
    dA0dY0(5,k+3) = + P_aM*dY0dA0(k  ,0) + P_eM*dY0dA0(k  ,1);
  };

  // State transition matrix
  dYdY0 = dYdA0 * dA0dY0;
    printf("\tdso::propagate_state Matrix dYdY0: (note dt=%.3f)\n",dt);
    for (int k=0;k<6;k++) {
      printf("\t\t|");
      for (int j=0; j<6; j++)
        printf("%8.1f ", dYdY0(k,j));
      printf("\n");
    }

  // propagate state
  return dso::propagate_state(GM, elements, dt);
}
