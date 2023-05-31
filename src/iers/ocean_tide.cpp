#include "tides.hpp"
#include <stdexcept>

dso::OceanTide::OceanTide(const char *fn, double scale, int max_degree,
                          int max_order, double GM, double R)
    : dCS(max_degree, max_order, GM, R) {
  /* parse file and collect main waves */
  if (dso::memmap_octide_coefficients(fn, doodsonFreqs, max_degree, max_order,
                                      scale)) {
    throw std::runtime_error(
        "ERROR Failed to parse Ocean Tide coefficients file\n");
  }
}
