#include "tides.hpp"
#include <stdexcept>

dso::AtmosphericTide::AtmosphericTide(const char *fn, double scale, int max_degree,
                          int max_order, double GM, double R)
    : dCS(max_degree, max_order, GM, R), V(max_degree + 3, max_degree + 3),
      W(max_degree + 3, max_degree + 3) {
  /* parse file and collect main waves */
  if (dso::memmap_octide_coefficients(fn, doodsonFreqs, max_degree, max_order,
                                      scale)) {
    throw std::runtime_error(
        "ERROR Failed to parse Atmospheric Tide coefficients file\n");
  }
}
