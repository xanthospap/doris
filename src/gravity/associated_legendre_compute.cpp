#include "associated_legendre.hpp"
#include <cmath>

namespace {
  // Associated Legendre Polynomials up to (n,m)=(4,4) according to Wikipedia
  // see https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#The_first_few_associated_Legendre_functions

void fill_degree4x4(
    double x,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &P) noexcept 
    {
      const double rpmx = std::sqrt(1e0 -x*x);
      const double pmx = rpmx * rpmx;
      const double x2 = x*x;

      P(0, 0) = 1e0;

      P(1, 0) = x;
      P(1, 1) = -rpmx;

      P(2, 0) = .5e0 * (3e0 * x2 - 1e0);
      P(2, 1) = -3e0 * x * rpmx;
      P(2, 2) = 3e0 * pmx;

      P(3, 0) = .5e0 * (5e0 * x * x2 - 3e0 * x);
      P(3, 1) = (3e0 / 2e0) * (1e0 - 5 * x2) * rpmx;
      P(3, 2) = 15e0 * x * pmx;
      P(3, 3) = -15e0 * rpmx * pmx;

      P(4, 0) = (1e0 / 8e0) * (35e0 * x2 * x2 - 30e0 * x2 + 3e0);
      P(4, 1) = -(5e0 / 2e0) * (7e0 * x2 * x - 3e0 * x) * rpmx;
      P(4, 2) = (15e0 / 2e0) * (7e0 * x2 - 1e0) * pmx;
      P(4, 3) = -105e0 * x * rpmx * pmx;
      P(4, 4) = 105e0 * pmx * pmx;

      return;
    }

}// unnamed namespace

/*
void dso::AssociatedLegendreFunctions::compute(double x) noexcept {
  // fill polynomials to min degree/order = (4,4)
  fill_degree4x4(x,P);

  for (int m = 1; m <= m_degree; m++) {
    P(m, m) = static_cast<double>(2 * m - 1) * cf * P(m - 1, m - 1);
    // P(m,m-1) = static_cast<double>(2*(m-1)+1) * sf * P(m,m);
  }
  for (int m = 0; m <= m_degree; m++) {
    P(m + 1, m) = static_cast<double>(2 * m + 1) * sf * P(m, m);
  }
  for (int n = 2; n <= m_degree; n++) {
    for (int m = 0; m < n + 1 - 2; m++) {
      P(n, m) = (1e0 / (n - m)) *
                ((2 * n - 1) * sf * P(n - 1, m) - (n + m - 1) * P(n - 2, m));
    }
  }

  return;
}
*/
