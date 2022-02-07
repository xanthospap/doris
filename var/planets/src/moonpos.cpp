#include "celgeo/iau.hpp"
#include "geodesy/units.hpp"
#include "planetpos.hpp"
#ifdef DEBUG
#include <cstdio>
#endif

using std::cos;
using std::sin;

int dso::moon_vector(double t, double *rsun) noexcept {

  // fundamental arguments using IAU conventions:
  const double Omega = iers2010::sofa::faom03(t);
  // mean anomaly of the Moon.
  const double l = iers2010::sofa::fal03(t);
  // mean anomaly of the Sun.
  const double lp = iers2010::sofa::falp03(t);
  // mean longitude of the Moon minus mean longitude of the ascending node
  const double F = iers2010::sofa::faf03(t);
  // mean elongation of the Moon from the Sun.
  const double D = iers2010::sofa::fad03(t);
  // mean longitude of the moon
  const double L0 = Omega + F;

#ifdef DEBUG
  printf("Fundamental Arguments for t=%15.5f\n", t);
  printf("L0 = %15.10e\n", L0);
  printf("l  = %15.10e\n", l);
  printf("lp = %15.10e\n", lp);
  printf("D  = %15.10e\n", D);
  printf("F  = %15.10e\n", F);
#endif

  const double c[] = {22640e0, 769e0, -4586e0, 2370e0, -668e0, -412e0, -212e0,
                      -206e0,  192e0, -165e0,  148e0,  -125e0, -110e0, -55e0};

  // lunar longitude
  const double lambda =
      L0 + c[0] * sin(l) + c[1] * sin(2e0 * l) + c[2] * sin(l - 2e0 * D) +
      c[3] * sin(2e0 * D) + c[4] * sin(lp) - c[5] * sin(2e0 * F) +
      c[6] * sin(2e0 * l - 2e0 * D) - c[7] * sin(l + lp - 2e0 * D) +
      c[8] * sin(l + 2e0 * D) - c[9] * sin(lp - 2e0 * D) + c[10] * sin(l - lp) -
      c[11] * sin(D) + c[12] * sin(l + lp) - c[13] * sin(2e0 * F - 2e0 * D);

  // lunar latitude
  const double beta =
      18520e0 * sin(F + lambda - L0 + 412e0 * sin(2e0 * F) + 541e0 * sin(lp)) -
      526e0 * sin(F - 2e0 * D) + 44e0 * sin(l + F - 2e0 * D) -
      31e0 * sin(-l + F - 2e0 * D) - 25e0 * sin(-2e0 * l + F) -
      23e0 * sin(lp + F - 2e0 * D) + 21e0 * sin(-l + F) +
      11e0 * sin(-lp + F - 2e0 * D);

  // distance
  const double r = 385000e0 - 20905e0 * cos(l) - 3699e0 * cos(2e0 * D - l) -
                   2956e0 * cos(2e0 * D) - 570e0 * cos(2e0 * l) +
                   246e0 * cos(2e0 * l - 2e0 * D) - 205e0 * cos(lp - 2e0 * D) -
                   171e0 * cos(l + 2e0 * D) - 152e0 * cos(l + lp - 2e0 * D);

  const double x = r * cos(lambda) * cos(beta);
  const double y = r * sin(lambda) * cos(beta);
  const double z = r * sin(beta);
  constexpr double ecliptic_obliquity = dso::deg2rad<double>(23.43929111);
  const double C = cos(ecliptic_obliquity);
  const double S = sin(ecliptic_obliquity);

  rsun[0] = x;
  rsun[1] = C * y + S * z;
  rsun[2] = -S * y + C * z;

  return 0;
}
