#include "associated_legendre.hpp"
#include <cmath>

void dso::AssociatedLegendreFunctions::compute(double phi) noexcept {
  const double sf = std::sin(phi);
  const double cf = std::cos(phi);
  P(0, 0) = 1e0;
  for (int m = 1; m <= m_order; m++) {
    P(m, m) = static_cast<double>(2 * m - 1) * cf * P(m - 1, m - 1);
    // P(m,m-1) = static_cast<double>(2*(m-1)+1) * sf * P(m,m);
  }
  for (int m = 0; m <= m_order; m++) {
    P(m + 1, m) = static_cast<double>(2 * m + 1) * sf * P(m, m);
  }
  for (int n = 2; n <= m_order; n++) {
    for (int m = 0; m < n + 1 - 2; m++) {
      P(n, m) = (1e0 / (n - m)) *
                ((2 * n - 1) * sf * P(n - 1, m) - (n + m - 1) * P(n - 2, m));
    }
  }
}