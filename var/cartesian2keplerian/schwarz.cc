#include "orbcrd.hpp"
#include <cstdio>

/// https://downloads.rene-schwarz.com/download/M002-Cartesian_State_Vectors_to_Keplerian_Orbit_Elements.pdf
void cartesian2keplerian_schwarz(const double *x, const double *v) noexcept {
  using math_3d::cross_product;
  using math_3d::dot_product;
  /*
   * Output:
   * -----------------------------------
   * semi-major axis, a (m)
   * eccentricity, e
   * argument of periapsis, ω (rad)
   * longitude of ascending node, Ω (rad)
   * inclination, i (rad)
   * mean anomaly, M (rad)
   */

  // Calculate orbital momentum vector
  double hv[3];
  cross_product(x, v, hv);
  const double hsq = dot_product(hv, hv);
  const double h = std::sqrt(hsq);
  const double rsq = dot_product(x, x);
  const double r = std::sqrt(rsq);

  // eccentricity vector
  double ev[3];
  cross_product(v, hv, ev);
  ev[0] /= ngpt::GM;
  ev[1] /= ngpt::GM;
  ev[2] /= ngpt::GM;
  ev[0] -= (x[0] / r);
  ev[1] -= (x[1] / r);
  ev[2] -= (x[2] / r);
  const double esq = dot_product(ev, ev);

  // eccentricity, e
  double e = std::sqrt(esq);

  // Determine the vector n pointing towards the ascending node
  double nv[] = {-hv[1], hv[0], 0e0};
  double n = std::sqrt(dot_product(nv, nv));

  //  True anomaly
  double ta = std::acos(dot_product(ev, x) / e / r);
  if (dot_product(x, v) < 0e0)
    ta = ngpt::D2PI - ta;

  // Orbital inclination, i
  double i = std::acos(hv[2] / h);

  // Eccentric anomaly, E
  double E =
      2e0 * std::atan2(std::tan(ta / 2e0), std::sqrt((1e0 + e) / (1e0 - e)));
  // if (E<0e0) E += ngpt::D2PI;

  // longtitude of the acending node, Ω
  double Omega = std::acos(nv[0] / n);
  // if (nv[1]<0e0) Omega = ngpt::D2PI - Omega;

  // argument of periapsis, ω
  double omega = std::acos(dot_product(nv, ev) / n / e);
  // if (ev[2]<0e0) omega = ngpt::D2PI - omega;

  // Mean anomaly, M
  double M = E - e * std::sin(E);

  // Semi-major axis, a
  double a = 1e0 / ((2e0 / r) - (dot_product(v, v) / ngpt::GM));

  /*
  printf("a    = %15.6f\n", a);
  printf("e    = %15.13f\n", e);
  printf("ω    = %15.12f\n", ngpt::rad2deg(omega));
  printf("Ω    = %15.12f\n", ngpt::rad2deg(Omega));
  printf("i    = %15.12f\n", ngpt::rad2deg(i));
  printf("M    = %15.12f\n", ngpt::rad2deg(M));
  */
  return;
}
