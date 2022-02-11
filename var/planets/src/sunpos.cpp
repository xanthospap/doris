#include "celgeo/iau.hpp"
#include "geodesy/units.hpp"
#include "planetpos.hpp"

using std::cos;
using std::sin;

// double Frac(double x) noexcept { return x - floor(x); }
constexpr double pi2 = iers2010::D2PI;

int dso::sun_vector_approx(double t, double *rsun) noexcept {
  double ipart; // dummy var

  // Mean anomaly [rad], ecliptic longitude [rad] and radius [m]
  const double M = pi2 * std::modf(0.9931267e0 + 99.9973583e0 * t, &ipart);
  const double L =
      pi2 * std::modf(0.7859444e0 + M / pi2 +
                          (6892e0 * sin(M) + 72e0 * sin(2e0 * M)) / 1296.0e3,
                      &ipart);
  const double r = 149.619e9 - 2.499e9 * cos(M) - 0.021e9 * cos(2 * M);

  // auxilary values [km]
  const double x = (r * cos(L)) / 1e3;
  const double y = (r * sin(L)) / 1e3;

  constexpr double ecliptic_obliquity = dso::deg2rad<double>(23.43929111e0);
  constexpr double C = cos(-ecliptic_obliquity);
  constexpr double S = sin(-ecliptic_obliquity);

  // X = R_x(-epsilon) * Vec
  rsun[0] = x;
  rsun[1] = C * y;
  rsun[2] = -S * y;

  return 0;
}
