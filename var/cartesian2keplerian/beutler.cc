#include "orbcrd.hpp"
#include <cstdio>

void cartesian2keplerian_beutler(const double *x, const double *xdot,
                                 double t) noexcept {
  using math_3d::dot_product;
  using math_3d::cross_product;
  /*
   * Output:
   * -----------------------------------
   * semi-latus rectum, p
   * eccentricity, e
   * argument of periapsis, ω (rad)
   * longitude of ascending node, Ω (rad)
   * inclination, i (rad)
   * time of perihelion, T0
   */
  double hv[3];
  cross_product(x, xdot, hv);
  hv[1] = -hv[1]; // bernese-like cross product
  double hsq = dot_product(hv, hv);
  double h = std::sqrt(hsq);

  // ascending node, Ω
  double Omega = std::atan2(hv[0], hv[1]);

  // inclination, i
  double i = std::acos(hv[2] / h);

  // constant Edot
  double rsq = dot_product(x, x);
  double r = std::sqrt(rsq);
  double Edot = dot_product(xdot, xdot) / 2e0 - ngpt::GM / r;

  // argument p
  double p = hsq / ngpt::GM;

  // argument e
  double e = std::sqrt(1e0 + 2e0 * hsq * Edot / ngpt::GM / ngpt::GM);

  // compute the argument of latitude at time t
  // this is: R_x(i) * R_z(Ω) * r (expanded to ...)
  double ck = std::cos(Omega);
  double sk = std::sin(Omega);
  double ci = std::cos(i);
  double si = std::sin(i);
  double xx0 = ck * x[0] + sk * x[1];
  double xx1 = -ci * sk * x[0] + ci * ck * x[1] + si * x[2];
  // double xx2 = si*sk*x[0]-si*ck*x[1]+ci*x[2];
  double u = std::atan2(xx1, xx0);
  double ecv = p / r - 1e0;
  // double esv = h / ngpt::GM;
  double esv = std::sqrt(p / ngpt::GM) / r * dot_product(x, xdot);
  double v1 = std::atan2(esv, ecv);
  double omega = u - v1;
  // for elliptic orbits (eq. 5.33)
  double a = p / (1e0 - e * e);
  double sinv = std::sin(v1);
  double cosv = std::cos(v1);
  /*
  double sinE = (std::sqrt(1e0-e*e) / (1e0+e)) * (sinv/cosv); // 4.51
  double cosE = (e + cosv) / (1e0+e*cosv);                    // 4.51
  double E = std::atan(sinE/cosE);
  */
  double E = std::atan(std::sqrt(1e0 + e * e) * sinv / (e + cosv));
  double T0 =
      t - (E - e * std::sin(E)) / std::sqrt(ngpt::GM / a / a / a); // Table 4.2

  /*
  printf("a    = %15.6f\n", a);
  printf("e    = %15.13f\n", e);
  printf("ω    = %15.12f\n", ngpt::rad2deg(omega));
  printf("Ω    = %15.12f\n", ngpt::rad2deg(Omega));
  printf("i    = %15.12f\n", ngpt::rad2deg(i));
  printf("T0   = %15.12f\n", T0);
  */
}
