#include "orbcrd.hpp"
#include <cstdio>

void cartesian2keplerian_bernese(const double *x, const double *v,
                                 double t) noexcept {
  using math_3d::dot_product;
  using math_3d::cross_product;
  /*
   * Output:
   * -----------------------------------
   * semi-major axis, a (m)
   * eccentricity, e
   * argument of periapsis, ω (rad) -- in bernese PER --
   * longitude of ascending node, Ω (rad) -- in bernese KN --
   * inclination, i (rad)
   * time of perigee/helion passing, T0
   */
  double h[3] = {
      x[1] * v[2] - x[2] * v[1],
      -x[2] * v[0] + x[0] * v[2],
      x[0] * v[1] - x[1] * v[0],
  };

  double kn = std::atan2(h[0], h[1]);

  double i = std::atan2(std::sqrt(h[0] * h[0] + h[1] * h[1]), h[2]);

  double p = (h[0] * h[0] + h[1] * h[1] + h[2] * h[2]) / ngpt::GM;
  double r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  double ecv = p / r - 1e0;
  double esv =
      std::sqrt(p / ngpt::GM) / r * (x[0] * v[0] + x[1] * v[1] + x[2] * v[2]);
  double v1 = std::atan2(esv, ecv);
  double e = std::sqrt(ecv * ecv + esv * esv);

  double ck = std::cos(kn);
  double sk = std::sin(kn);
  double ci = std::cos(i);
  double si = std::sin(i);
  double xx[] = {ck * x[0] + sk * x[1],
                 -ci * sk * x[0] + ci * ck * x[1] + si * x[2],
                 si * sk * x[0] - si * ck * x[1] + ci * x[2]};

  double u = std::atan2(xx[1], xx[0]);
  double per = u - v1;
  double ex = 2e0 * std::atan(std::sqrt((1e0 - e) / (1e0 + e)) *
                              (std::sin(v1 / 2) / std::cos(v1 / 2e0)));
  double a = p / (1e0 - e * e);
  double a3 = a * a * a;
  double t0 = t - (ex - e * std::sin(ex)) / std::sqrt(ngpt::GM / a3);

  /*
  printf("a    = %15.6f\n", a);
  printf("e    = %15.13f\n", e);
  printf("ω    = %15.12f\n", ngpt::rad2deg(per));
  printf("Ω    = %15.12f\n", ngpt::rad2deg(kn));
  printf("i    = %15.12f\n", ngpt::rad2deg(i));
  printf("T0   = %15.12f\n", t0);
  */
  return;
}
