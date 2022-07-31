#include "nrlmsise00.hpp"

/// @brief chemistry/dissociation correction for msis model
/// @param[in] alt Altitude
/// @param[in] r   Target ratio
/// @param[in] h1  Transition scale length
/// @param[in] zh  Altitude of 1/2 r
double dso::nrlmsise00::detail::ccor(double alt, double r, double h1,
                             double zh) noexcept {
  constexpr const double echeck = 70e0;
  const double e = (alt - zh) / h1;

  double ccor;
  if (e > echeck) {
    ccor = 0e0;
  } else if (e < -echeck) {
    ccor = r;
  } else {
    const double ex = std::exp(e);
    ccor = r / (1e0 + ex);
  }

  return std::exp(ccor);
}
