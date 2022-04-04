#include "astrodynamics.hpp"
#include <cmath>

double dso::kepler(double e, double M, double tolerance_rad) noexcept {
  double E = (e > 0.8e0) ? M_PI : M;
  double En;
  do {
    En = E - (E - e * std::sin(E) - M) / (1e0 - e * std::cos(E));
    E = En;
  } while (E - e * std::sin(E) - M > tolerance_rad);

  return E;
}
