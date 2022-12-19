#include "astrodynamics.hpp"
#include "iers2010/iersc.hpp"
#include <cmath>

constexpr const int max_it = 20;
constexpr const double pi = iers2010::DPI;
constexpr const double pi2 = iers2010::D2PI;

/// @brief Normalize in range [0, 2π], including negative angles
inline double arg02pi(double angle) noexcept {
  const double t = std::fmod(angle, pi2);
  return t + (t < 0e0) * pi2;
}

double dso::kepler(double e, double M, int &status,
                   double tolerance_rad) noexcept {
  status = 0;
  M = arg02pi(M);
  const bool elt08 = (e < 8e-1);
  double E = (!elt08) * pi + elt08 * M;

  double f;
  int it = 0;
  do {
    f = E - e * std::sin(E) - M;
    E = E - f / (1e0 - e * std::cos(E));
  } while (std::abs(f) > tolerance_rad && ++it < max_it);

#ifdef COUNT_KEPLER_ITERATIONS
  status = (it >= max_it) ? -1 : it;
#else
  status = it >= max_it;
#endif
  return E;
}

double dso::kepler_vallado(double e, double M, int &status,
                           double tolerance_rad) noexcept {
  status = 0;
  M = arg02pi(M);
  const bool mr = (M > -pi && M < 0) || (M > pi);
  const double esign = mr * -e + (!mr) * e;
  double En = M + esign;

  double E;
  int it = 0;
  do {
    E = En;
    En = E + (M - E + e * std::sin(E)) / (1e0 - e * std::cos(E));
  } while (std::abs(E - En) > tolerance_rad && ++it < max_it);

#ifdef COUNT_KEPLER_ITERATIONS
  status = (it >= max_it) ? -1 : it;
#else
  status = it >= max_it;
#endif
  return En;
}
