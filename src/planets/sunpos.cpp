#include "gcem.hpp" // for constexpr math (trigonometric funcs)
#include "geodesy/units.hpp"
#include "iers2010/iau.hpp"
#include "planetpos.hpp"
#ifdef DEBUG
#include <cstdio>
#endif

// Astronomical unit (km, IAU 2012)
#define DAUKM (149597870.7e0)

/// @brief Normalize a given angle (in degrees) to the range [0,1]
constexpr double normalize_deg(double angle) noexcept { return angle / 360e0; }
constexpr double pi2 = iers2010::D2PI;

using std::cos;
using std::sin;

/// Note that this implementation uses normalized angles (and not the
/// nominal values recorded in the book.
int dso::sun_vector_montenbruck(double t, double *rsun) noexcept {
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
  constexpr double C = gcem::cos(-ecliptic_obliquity);
  constexpr double S = gcem::sin(-ecliptic_obliquity);

  // X = R_x(-epsilon) * Vec
  rsun[0] = x;
  rsun[1] = C * y;
  rsun[2] = -S * y;

  return 0;
}

int dso::sun_vector_vallado(double t, double *rsun) noexcept {
  const double lambda =
      std::fmod(280.460e0 + 36000.771e0 * t, 360e0); // [degrees]
  const double M =
      std::fmod(357.5291092e0 + 35999.50534e0 * t, 360e0); // [degrees]
  const double Mrad = dso::deg2rad(M);
  const double l_ecliptic = std::fmod(lambda + 1.914666471e0 * sin(Mrad) +
                                          0.019994643e0 * sin(2e0 * Mrad),
                                      360e0);           // [degrees]
  const double epsilon = 23.439291e0 - 0.0130042e0 * t; // [degrees]
  const double r_au = 1.000140612e0 - 0.016708617e0 * cos(Mrad) -
                      0.000139589e0 * cos(2e0 * Mrad); // [degrees]

  const double S = sin(dso::deg2rad<double>(l_ecliptic));
  rsun[0] = r_au * cos(dso::deg2rad<double>(l_ecliptic));
  rsun[1] = r_au * cos(dso::deg2rad<double>(epsilon)) * S;
  rsun[2] = r_au * sin(dso::deg2rad<double>(epsilon)) * S;

  // AU to km
  for (int i = 0; i < 3; i++)
    rsun[i] *= DAUKM;

  return 0;
}

#ifdef DEBUG
// This is again the implementation (2), following Vallado, but with normalized
// angles; makes no difference in accuracy whatsoever. Have not tested speed
int dso::sun_vector_approx21(double t, double *rsun) noexcept {
  double dummy;
  constexpr double cf1 = normalize_deg(280.460e0);
  constexpr double cf2 = normalize_deg(36000.771e0);
  const double lambda = std::modf(cf1 + cf2 * t, &dummy) * pi2; // [rad]

  constexpr double cf3 = normalize_deg(357.5291092e0);
  constexpr double cf4 = normalize_deg(35999.50534e0);
  const double M = std::modf(cf3 + cf4 * t, &dummy) * pi2; // [rad]

  const double lon_arg = 1.914666471e0 * sin(M) + 0.019994643e0 * sin(2e0 * M);
  const double lon =
      std::fmod(normalize_deg(lon_arg) * pi2 + lambda, iers2010::D2PI);

  const double e = dso::deg2rad<double>(23.439291e0 - 0.0130042e0 * t); // [rad]
  const double r_au =
      1.000140612e0 - 0.016708617e0 * cos(M) - 0.000139589e0 * cos(2e0 * M);

  const double S = sin(lon);
  rsun[0] = r_au * cos(lon);
  rsun[1] = r_au * cos(e) * S;
  rsun[2] = r_au * sin(e) * S;

  // AU to km
  for (int i = 0; i < 3; i++)
    rsun[i] *= DAUKM;

  return 0;
}
#endif

int dso::sun_vector_curtis(double n, double *rsun) noexcept {
  // Calculate the mean anomaly M [deg]
  const double M = std::fmod(357.529e0 + 0.98560023e0 * n, 360e0);
  // Calculate the mean solar longitude L [deg]
  const double L = std::fmod(280.459e0 + 0.98564736e0 * n, 360e0);
  // Calculate the longitude λ [deg]
  const double lon = std::fmod(L + 1.915e0 * sin(dso::deg2rad(M)) +
                                   0.0200e0 * sin(2e0 * dso::deg2rad(M)),
                               360e0);
  // Calculate the obliquity ε [deg]
  const double e = 23.439e0 - 3.65e-7 * n;
  // Calculate the unit vector u_n from the earth to the sun
  const double ux = cos(dso::deg2rad(lon));
  const double uy = cos(dso::deg2rad(e)) * sin(dso::deg2rad(lon));
  const double uz = sin(dso::deg2rad(e)) * sin(dso::deg2rad(lon));
  // Calculate the distance r_S of the sun from the earth [AU]
  const double rs = 1.00014e0 - 0.01671e0 * cos(dso::deg2rad(M)) -
                    0.000140e0 * cos(2e0 * dso::deg2rad(M));
  // Calculate the sun’s geocentric position vector R_s = r_S * u_n
  rsun[0] = (ux * rs) * DAUKM;
  rsun[1] = (uy * rs) * DAUKM;
  rsun[2] = (uz * rs) * DAUKM;

  return 0;
}
