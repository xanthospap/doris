#include "nrlmsise00.hpp"

double dso::nrlmsise00::splint(const double *__restrict__ xa,
                               const double *__restrict__ ya,
                               double *__restrict__ y2a, int n,
                               double x) noexcept {

  int klo = 0, khi = n - 1;
  do {
    int k = (khi + klo) / 2;

    if (xa[k] > x) {
      khi = k;
    } else {
      klo = k;
    }

  } while (khi - klo > 1);

  const double h = xa[khi] - xa[klo];
  if (h == 0e0)
    return 1e-30;

  const double a = (xa[khi] - x) / h;
  const double b = (x - xa[klo]) / h;
  const double y =
      a * ya[klo] + b * ya[khi] +
      ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * h * h / 6e0;

  return y;
}
