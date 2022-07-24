#include "nrlmsise00.hpp"

/// @brief Chemistry/dissociation correction for msis models
/// @param[in] alt  Altitude
/// @param[in] r    Target ratio
/// @param[in] h1   Transition scale length
/// @param[in] zh   Altitude of 1 / 2 r
double dso::nrlmsise00::ccor2(double alt, double r, double h1, double zh,
                              double h2) noexcept {
  constexpr const double echeck = 70e0;

  const double e1 = (alt - zh) / h1;
  const double e2 = (alt - zh) / h2;

  double ccor2;
  if (e1 > echeck || e2 > echeck) {
    ccor2 = 0e0;
  } else if (e1 < -echeck && e2 < -echeck) {
    ccor2 = r;
  } else {
    const double ex1 = std::exp(e1);
    const double ex2 = std::exp(e2);
    ccor2 = r / (1e0 + 0.5e0 * (ex1 + ex2));
  }

  return std::exp(ccor2);
}
