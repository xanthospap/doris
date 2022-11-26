#include "egravity.hpp"
#include <cmath>
#include <cstdio>

int dso::lagrange_polynomials(
    double x, double y, double z, double R, int l, int k,
    dso::Mat2D<dso::MatrixStorageType::Trapezoid> &V,
    dso::Mat2D<dso::MatrixStorageType::Trapezoid> &W) noexcept {

  if (k > l) {
    fprintf(stderr,
            "[ERROR] Invalid degree/order given to compute Lagrange "
            "polynomials! given degre=%d order=%d (traceback: %s)\n",
            l, k, __func__);
    return 1;
  }

  // distance squared
  const double r2 = x * x + y * y + z * z;
  // rho factor, R^2 / r^2
  const double rho = R * R / r2;
  // Normalized coordinates
  const double x0 = R * x / r2;
  const double y0 = R * y / r2;
  const double z0 = R * z / r2;

  // Calculate zonal terms V(n,0); set W(n,0)=0.0
  V(0, 0) = R / std::sqrt(r2);
  W(0, 0) = 0e0;

  V(1, 0) = z0 * V(0, 0);
  W(1, 0) = 0e0;

  int n = 0, m = 0;
  for (n = 2; n <= l; n++) {
    V(n, 0) =
        ((2e0 * n - 1e0) * z0 * V(n - 1, 0) - (n - 1) * rho * V(n - 2, 0)) / n;
    W(n, 0) = 0e0;
  }

  for (m = 1; m <= k; m++) {
    // use 3.29 to compute V_mm and W_mm
    V(m, m) = (2 * m - 1) * (x0 * V(m - 1, m - 1) - y0 * W(m - 1, m - 1));
    W(m, m) = (2 * m - 1) * (x0 * W(m - 1, m - 1) + y0 * V(m - 1, m - 1));

    if (m < l) {
      V(m + 1, m) = (2 * m + 1) * z0 * V(m, m);
      W(m + 1, m) = (2 * m + 1) * z0 * W(m, m);
    };

    // use 3.30 to compute V_nm and W_nm
    for (n = m + 2; n <= l; n++) {
      V(n, m) =
          ((2 * n - 1) * z0 * V(n - 1, m) - (n + m - 1) * rho * V(n - 2, m)) /
          (n - m);
      W(n, m) =
          ((2 * n - 1) * z0 * W(n - 1, m) - (n + m - 1) * rho * W(n - 2, m)) /
          (n - m);
    }
  }

  return 0;
}
