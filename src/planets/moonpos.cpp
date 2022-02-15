#include "iers2010/iau.hpp"
#include "geodesy/units.hpp"
#include "planetpos.hpp"
#ifdef DEBUG
#include <cstdio>
#endif

using std::cos;
using std::sin;

constexpr double pi2 = iers2010::D2PI;

int dso::moon_vector_approx(double t, double *rsun) noexcept {
  double ipart; // dummy var

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
  const double L0 = pi2 * std::modf(0.606433e0 + 1'336.851344e0 * t, &ipart);

  // lunar longitude correction [sec]
  const double dlambda =
      22'640e0 * sin(l) + 769e0 * sin(2e0 * l) - 4'586e0 * sin(l - 2e0 * D) +
      2'370e0 * sin(2e0 * D) + -668e0 * sin(lp) - 412e0 * sin(2e0 * F) +
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
      std::fmod(18'520e0 * sin(sinarg) - 526e0 * sin(F - 2e0 * D) +
                    44e0 * sin(l + F - 2e0 * D) - 31e0 * sin(-l + F - 2e0 * D) -
                    25e0 * sin(-2e0 * l + F) - 23e0 * sin(lp + F - 2e0 * D) +
                    21e0 * sin(-l + F) + 11e0 * sin(-lp + F - 2e0 * D),
                iers2010::TURNAS) *
      iers2010::DAS2R;

  // distance [m]
  const double r = (385'000e0 - 20'905e0 * cos(l) - 3'699e0 * cos(2 * D - l) -
                    2'956e0 * cos(2 * D) - 570e0 * cos(2 * l) +
                    246e0 * cos(2 * l - 2 * D) - 205e0 * cos(lp - 2 * D) -
                    171e0 * cos(l + 2 * D) - 152e0 * cos(l + lp - 2 * D)) *
                   1e3;

  // auxilary values [km]
  const double x = (r * cos(lambda) * cos(beta)) / 1e3;
  const double y = (r * sin(lambda) * cos(beta)) / 1e3;
  const double z = (r * sin(beta)) / 1e3;

  constexpr double ecliptic_obliquity = dso::deg2rad<double>(23.43929111e0);
  const double C = cos(-ecliptic_obliquity);
  const double S = sin(-ecliptic_obliquity);

  // X = R_x(-epsilon) * Vec
  rsun[0] = x;
  rsun[1] = C * y + S * z;
  rsun[2] = -S * y + C * z;

  return 0;
}

int dso::moon_vector_vallado(double t, double *rmoon) noexcept {
  const double lecliptic =
      218.32e0 + 481'267.8813e0 * t +
      6.29e0 * sin(dso::deg2rad<double>(134.9e0 + 477'198.85e0 * t)) -
      1.27e0 * sin(dso::deg2rad<double>(259.2e0 - 413'335.38e0 * t)) +
      0.66e0 * sin(dso::deg2rad<double>(235.7e0 + 890'534.23e0 * t)) +
      0.21e0 * sin(dso::deg2rad<double>(269.9e0 + 954'397.70 * t)) -
      0.19e0 * sin(dso::deg2rad<double>(357.5e0 + 35'999.05e0 * t)) -
      0.11e0 * sin(dso::deg2rad<double>(186.6e0 + 966'404.05e0 * t)); // [deg]

  const double fecliptic =
      5.13e0 * sin(dso::deg2rad<double>(93.3e0 + 483'202.03e0 * t)) +
      0.28e0 * sin(dso::deg2rad<double>(228.2e0 + 960'400.87e0 * t)) -
      0.28e0 * sin(dso::deg2rad<double>(318.3e0 + 6'003.18e0 * t)) -
      0.17e0 * sin(dso::deg2rad<double>(217.6e0 - 407'332.20e0 * t)); // [deg]

  const double rho =
      0.9508e0 + 0.0518e0 * cos(dso::deg2rad<double>(134.9e0 + 477'198.85 * t)) +
      0.0095e0 * cos(dso::deg2rad<double>(259.2e0 - 413'335.38e0 * t)) +
      0.0078e0 * cos(dso::deg2rad<double>(235.7e0 + 890'534.23e0 * t)) +
      0.0028e0 * cos(dso::deg2rad<double>(269.9e0 + 954'397.70e0 * t)); // [deg]

  const double obliquity =
      23.439291e0 + t * (-0.0130042e0 + t * (-1.64e-7 + 5.04e-7 * t)); // [deg]

  // transform to radians ....
  const double L = std::fmod(dso::deg2rad(lecliptic), iers2010::D2PI);
  const double F = std::fmod(dso::deg2rad(fecliptic), iers2010::D2PI);
  const double R = std::fmod(dso::deg2rad(rho), iers2010::D2PI);
  const double E = std::fmod(dso::deg2rad(obliquity), iers2010::D2PI);
  // compute triginometric numbers ....
  const double Cf = cos(F);
  const double Sf = sin(F);
  const double Sl = sin(L);
  const double Cl = cos(L);
  const double Ce = cos(E);
  const double Se = sin(E);
  // geocentric direction cosines ...
  const double l = Cf * Cl;
  const double m = Ce * Cf * Sl - Se * Sf;
  const double n = Se * Cf * Sl + Ce * Sf;
  // moon position vector in km
  const double Rm = /*iers2010::Re*/6378136.6e0 / 1e3 / sin(R); // [km]
  rmoon[0] = Rm * l;
  rmoon[1] = Rm * m;
  rmoon[2] = Rm * n;

  return 0;
}
