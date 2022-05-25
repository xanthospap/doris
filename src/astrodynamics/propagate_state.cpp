#include "astrodynamics.hpp"
#include <iers2010/matvec.hpp>
#include <cstdio>

using dso::Vector3;

int dso::propagate_state(double GM, const Vector3 &r0, const Vector3 &v0,
                         double dt, Vector3 &r, Vector3 &v) noexcept {
  
  // state to orbital elements
  dso::OrbitalElements ele;
  int error = state2elements(r0, v0, ele, GM);
  //printf("\tcomputed elements: ");
  //printf("%.5f ", ele.semimajor());
  //printf("%.5f ", ele.eccentricity());
  //printf("%.5f ", ele.inclination());
  //printf("%.5f ", ele.Omega());
  //printf("%.5f ", ele.omega());
  //printf("%.5f ", ele.mean_anomaly());
  //printf("\n");

  // propagate state to next epoch (dt)
  error = dso::elements2state(ele, dt, r, v, GM);

  return error;
}
