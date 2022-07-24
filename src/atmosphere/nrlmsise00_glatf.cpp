#include "nrlmsise00.hpp"
#include "geodesy/units.hpp"

using namespace dso::nrlmsise00;

/// @brief calculate latitude variable gravity (gv) and effective radius (reff)
double dso::Nrlmsise00::glatf(double lat, double &gv) const noexcept {
  const double c2 = std::cos(2e0 * dso::deg2rad(lat));
  gv = egrav * 1e2 * (1e0 - 0.0026373e0 * c2);
  return 2e0 * gsurf / (3.085462e-6 + 2.27e-9 * c2) * 1e-5;
}
