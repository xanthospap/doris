#include "nrlmsise00.hpp"
#include <algorithm>

double dso::nrlmsise00::splini(const double *__restrict__ xa,
                               const double *__restrict__ ya,
                               double *__restrict__ y2a, int n,
                               double x) noexcept {

  int klo = 0, khi = 1;
  double xx, a, a2, b, b2, h, yi;

  do {
    xx = x;
    if (khi < n - 1)
      xx = std::min(x, xa[khi]);
    h = xa[khi] - xa[klo];
    a = (xa[khi] - xx) / h;
    b = (xx - xa[klo]) / h;
    a2 = a * a;
    b2 = b * b;
    yi = yi + ((1e0 - a2) * ya[klo] / 2e0 + b2 * ya[khi] / 2e0 +
               ((-(1e0 + a2 * a2) / 4e0 + a2 / 2e0) * y2a[klo] +
                (b2 * b2 / 4e0 - b2 / 2e0) * y2a[khi]) *
                   h * h / 6e0) *
                  h;
    ++klo;
    ++khi;
  } while (x > xa[klo] && khi < n);

  return yi;
}
