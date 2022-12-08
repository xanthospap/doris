#include "cmat2d.hpp"
#include "egravity.hpp"
#include "geodesy/geodesy.hpp"
#include "harmonic_coeffs.hpp"
#include <array>
#include <cmath>
#ifdef DEBUG
#include <cassert>
#include <cfenv>
#include <limits>
//#  pragma STDC_FENV_ACCESS on
#endif

class Factor {
public:
  double f1(int n, int m) const noexcept {
  #ifdef DEBUG
    assert(n>0 && m<=n);
  #endif
    if (n == m)
      return (n == 1) ? std::sqrt(3.0) : std::sqrt((2. * n + 1.) / (2. * n));
    else
      return std::sqrt((2. * n + 1.) / static_cast<double>((n + m) * (n - m)) *
                (2. * n - 1.));
  }
  double f2(int n, int m) const noexcept {
    return -std::sqrt((2. * n + 1.) / static_cast<double>((n + m) * (n - m)) *
                      (n - m - 1.) * (n + m - 1.) / (2. * n - 3.));
  }
};

int test::gravacc0(const dso::HarmonicCoeffs &cs,
                   const Eigen::Matrix<double, 3, 1> &p, int degree, double Re,
                   double GM, Eigen::Matrix<double, 3, 1> &acc) noexcept {

  const int lp_degree = degree + 1; // aka, [0,....degree+1]
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> W(lp_degree + 2,
                                                            lp_degree + 2);
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> M(lp_degree + 2,
                                                            lp_degree + 2);
  W.fill_with(0e0);
  M.fill_with(0e0);

  const Factor F;

  {
    const Eigen::Matrix<double, 3, 1> r = p / Re;
    const double rr = std::pow(1 / r.norm(), 2);
    const double x = r.x() * rr;
    const double y = r.y() * rr;
    const double z = r.z() * rr;
    M(0, 0) = 1e280 / r.norm();

    // Recursion diagonal
    // C(n-1,n-1) -> C(n,n)
    for (int n = 1; n <= lp_degree; n++) {
      M(n, n) = F.f1(n, n) * (x * M(n - 1, n - 1) - y * W(n - 1, n - 1));
      W(n, n) = F.f1(n, n) * (y * M(n - 1, n - 1) + x * W(n - 1, n - 1));
    }

    // Recursion secondary diagonal
    // C(n-1,n-1) -> C(n,n-1)
    for (int n = 1; n <= lp_degree; n++) {
      M(n, n - 1) = F.f1(n, n - 1) * z * M(n - 1, n - 1);
      W(n, n - 1) = F.f1(n, n - 1) * z * W(n - 1, n - 1);
    }

    // Recursion others
    // C(n-1,m),C(n-1,m) -> C(n,m)
    for (int m = 0; m <= lp_degree; m++) {
      for (int n = m + 2; n <= lp_degree; n++) {
        M(n, m) = F.f1(n, m) * z * M(n - 1, m) + F.f2(n, m) * rr * M(n - 2, m);
        W(n, m) = F.f1(n, m) * z * W(n - 1, m) + F.f2(n, m) * rr * W(n - 2, m);
      }
    }

    M.multiply(1e-280);
    W.multiply(1e-280);
  }

  // acceleration in cartesian components
  acc = Eigen::Matrix<double, 3, 1>::Zero();
  const int minDegree = 1;

  for (int n = minDegree; n <= degree; n++) {
    // 0. Order
    double wm0 = std::sqrt(static_cast<double>(n + 1) * (n + 1));
    double wp1 =
        std::sqrt(static_cast<double>(n + 1) * (n + 2)) / std::sqrt(2.0);

    double Cm0 = wm0 * M(n + 1, 0);
    double Cp1 = wp1 * M(n + 1, 1);
    double Sp1 = wp1 * W(n + 1, 1);

    double gx = cs.C(n, 0) * (-2 * Cp1);
    double gy = cs.C(n, 0) * (-2 * Sp1);
    double gz = cs.C(n, 0) * (-2 * Cm0);
    
    for (int m = 1; m <= n; m++) {
      double wm1 = std::sqrt(static_cast<double>(n - m + 1) * (n - m + 2)) *
                   ((m == 1) ? std::sqrt(2.0) : 1.0);
      wm0 = std::sqrt(static_cast<double>(n - m + 1) * (n + m + 1));
      wp1 = std::sqrt(static_cast<double>(n + m + 1) * (n + m + 2));

      const double Cm1 = wm1 * M(n + 1, m - 1);
      const double Sm1 = wm1 * W(n + 1, m - 1);
      Cm0 = wm0 * M(n + 1, m);
      const double Sm0 = wm0 * W(n + 1, m);
      Cp1 = wp1 * M(n + 1, m + 1);
      Sp1 = wp1 * W(n + 1, m + 1);

      gx += cs.C(n, m) * (Cm1 - Cp1) + cs.S(n, m) * (Sm1 - Sp1);
      gy += cs.C(n, m) * (-Sm1 - Sp1) + cs.S(n, m) * (Cm1 + Cp1);
      gz += cs.C(n, m) * (-2 * Cm0) + cs.S(n, m) * (-2 * Sm0);
    }

    acc += std::sqrt((2. * n + 1.) / (2. * n + 3.)) *
           Eigen::Matrix<double, 3, 1>(gx, gy, gz);
  }

  acc *= GM / (2e0 * Re * Re);

  return 0;
}
