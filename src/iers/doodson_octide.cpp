#include "tides.hpp"
#include <stdexcept>

dso::DoodsonOceanTideConstituent::DoodsonOceanTideConstituent(
    const dso::DoodsonNumber d, int max_degree, int max_order)
    : doodson(d), maxl(max_degree), maxm(max_order),
      DelCpl(max_degree + 1, max_degree + 1),
      DelSpl(max_degree + 1, max_degree + 1),
      DelCmi(max_degree + 1, max_degree + 1),
      DelSmi(max_degree + 1, max_degree + 1) {
  if (max_order > max_degree || max_degree <= 1)
    throw std::runtime_error("[ERROR] Invalid degree/order sizes for "
                             "dso::DoodsonOceanTideConstituent\n");
}

void dso::DoodsonOceanTideConstituent::clear_coefficients() noexcept {
    DelCpl.fill_with(0e0);
    DelSpl.fill_with(0e0);
    DelCmi.fill_with(0e0);
    DelSmi.fill_with(0e0);
}

void dso::DoodsonOceanTideConstituent::resize(int maxDegree) noexcept {
    assert(maxDegree > 0);
    DelCpl.resize(maxDegree+1, maxDegree+1);
    DelSpl.resize(maxDegree+1, maxDegree+1);
    DelCmi.resize(maxDegree+1, maxDegree+1);
    DelSmi.resize(maxDegree+1, maxDegree+1);
    maxl = maxDegree;
    maxm = maxDegree;
}
