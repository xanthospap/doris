#include "celgeo/iau.hpp"
#include "geodesy/units.hpp"
#include "planetpos.hpp"
#ifdef DEBUG
#include <cstdio>
#endif

using std::cos;
using std::sin;

double Frac(double x) noexcept { return x - floor(x); }
constexpr double pi2 = iers2010::D2PI;

int dso::moon_vector(double t, double *rsun) noexcept {

  // fundamental arguments using IAU conventions [rad]:
  // mean anomaly of the Moon.
  const double l = iers2010::sofa::fal03(t);
  // mean anomaly of the Sun.
  const double lp = iers2010::sofa::falp03(t);
  // mean longitude of the Moon minus mean longitude of the ascending node
  const double F = iers2010::sofa::faf03(t);
  // mean elongation of the Moon from the Sun.
  const double D = iers2010::sofa::fad03(t);
  // mean longitude of the moon w.r.t. J2000 equinox
  const double L0 = pi2 * Frac(0.606433 + 1336.851344 * t);

  // lunar longitude correction [sec]
  const double dlambda =
      22640e0 * sin(l) + 769e0 * sin(2e0 * l) - 4586e0 * sin(l - 2e0 * D) +
      2370e0 * sin(2e0 * D) + -668e0 * sin(lp) - 412e0 * sin(2e0 * F) +
      -212e0 * sin(2e0 * l - 2e0 * D) - 206e0 * sin(l + lp - 2e0 * D) +
      192e0 * sin(l + 2e0 * D) - 165e0 * sin(lp - 2e0 * D) +
      148e0 * sin(l - lp) - 125e0 * sin(D) - 110e0 * sin(l + lp) -
      55e0 * sin(2e0 * F - 2e0 * D);

  // lunar longitude [rad]
  const double lambda =
      std::fmod(L0 * iers2010::DR2AS + dlambda, iers2010::TURNAS) *
      iers2010::DAS2R;

  // lunar latitude in [rad]
  const double sinarg = std::fmod(
      F + (dlambda + 412e0 * sin(2e0 * F) + 541e0 * sin(lp)) * iers2010::DAS2R,
      iers2010::D2PI);
  const double beta =
      std::fmod(18520e0 * sin(sinarg) - 526e0 * sin(F - 2e0 * D) +
                    44e0 * sin(l + F - 2e0 * D) - 31e0 * sin(-l + F - 2e0 * D) -
                    25e0 * sin(-2e0 * l + F) - 23e0 * sin(lp + F - 2e0 * D) +
                    21e0 * sin(-l + F) + 11e0 * sin(-lp + F - 2e0 * D),
                iers2010::TURNAS) *
      iers2010::DAS2R;

  // distance [m]
  const double r = (385000e0 - 20905e0 * cos(l) - 3699e0 * cos(2 * D - l) -
                    2956e0 * cos(2 * D) - 570e0 * cos(2 * l) +
                    246e0 * cos(2 * l - 2 * D) - 205e0 * cos(lp - 2 * D) -
                    171e0 * cos(l + 2 * D) - 152e0 * cos(l + lp - 2 * D)) *
                   1e3;

  // auxilary values
  const double x = r * cos(lambda) * cos(beta);
  const double y = r * sin(lambda) * cos(beta);
  const double z = r * sin(beta);

  constexpr double ecliptic_obliquity = dso::deg2rad<double>(23.43929111e0);
  const double C = cos(-ecliptic_obliquity);
  const double S = sin(-ecliptic_obliquity);

  // X = R_x(-epsilon) * Vec
  rsun[0] = x;
  rsun[1] = C * y + S * z;
  rsun[2] = -S * y + C * z;

  return 0;
}
