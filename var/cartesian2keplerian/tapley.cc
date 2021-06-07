#include "orbcrd.hpp"
#include <cstdio>

void cartesian2keplerian_tapley(const double *x, const double *xdot) noexcept {
  using math_3d::dot_product;
  using math_3d::cross_product;

  const double rsq = dot_product(x, x);
  const double r = std::sqrt(rsq);
  double hv[3];
  cross_product(x, xdot, hv);
  const double hsq = dot_product(hv, hv);
  const double h = std::sqrt(hsq);

  double i = std::acos(hv[2] / h);
  double Omega = std::atan2(hv[0], -hv[1]);

  double ksi = dot_product(xdot, xdot)/2e0 - ngpt::GM / r;
  double a = -ngpt::GM / ksi / 2e0;
  double e = std::sqrt(1e0 + (2e0*ksi*hsq) / ngpt::GM / ngpt::GM);
  double p = hsq / ngpt::GM;

  double cosf = (p-r) / r / e;
  double sinf = (p/h/e) * dot_product(x, xdot) / r;
  double cosof = (x[0]/r)*std::cos(Omega) + (x[1]/r)*std::sin(Omega);
  double sinof = x[2] / r / std::sin(i);
  double fangle = std::atan(sinf / cosf);
  double ofangle = std::atan(sinof / cosof);
  double omega = ofangle - fangle;

  double cosE = r * cosf / a + e;
  double b = a * std::sqrt(1e0 - e*e);
  double sinE = r * sinf / b;
  double E = std::atan(sinE / cosE);
  double M = E - e * sinE;
  
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
