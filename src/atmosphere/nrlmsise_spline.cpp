#include "nrlmsise00.hpp"

/// @brief Calculate 2nd derivatives of cubic spline interp function
void dso::nrlmsise00::spline(const double *__restrict__ x,
                             const double *__restrict__ y, int n, double yp1,
                             double ypn, double *__restrict__ y2,
                             double *work /*size >= n */) noexcept {
  double *__restrict__ u = work;

  if (yp1 > 0.99e30) {
    y2[0] = 0e0;
    u[0] = 0e0;
  } else {
    y2[0] = -0.5e0;
    u[0] = (3e0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  }

  for (int i = 1; i < n - 1; i++) {
    const double sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
    const double p = sig * y2[i - 1] + 2e0;
    y2[i] = (sig - 1e0) / p;
    u[i] = (6e0 *
                ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
                 (y[i] - y[i - 1]) / (x[i] - x[i - 1])) /
                (x[i + 1] - x[i - 1]) -
            sig * u[i - 1]) /
           p;
  }

  double qn, un;
  if (ypn > 0.99e30) {
    qn = un = 0e0;
  } else {
    qn = 0.5e0;
    const int nn = n - 1;
    un = (3e0 / (x[nn] - x[nn - 1])) *
         (ypn - (y[nn] - y[nn - 1]) / (x[nn] - x[nn - 1]));
  }

  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1e0);
  for (int k = n - 2; k >= 0; k--) {
    y2[k] = y2[k] * y2[k + 1] + u[k];
  }

  return;
}
