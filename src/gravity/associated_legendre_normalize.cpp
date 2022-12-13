#include "associated_legendre.hpp"
#include <cmath>

void dso::AssociatedLegendreFunctions::normalize() noexcept {
  const double sqrt2_ = std::sqrt(2e0);
  for (int n = 1; n <= m_degree; n++) {
    const double tnp1 = 2 * n + 1;
    double factor = std::sqrt(tnp1);
    P(n, 0) *= factor;
    factor *= sqrt2_;
    for (int m = 1; m <= n; m++) {
      factor /= std::sqrt(static_cast<double>(n + m));
      factor /= std::sqrt(static_cast<double>(n - (m - 1)));
      P(n, m) *= factor;
    }
  }
  return;
}
