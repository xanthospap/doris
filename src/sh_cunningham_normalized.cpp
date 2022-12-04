#include "cmat2d.hpp"
#include "harmonic_coeffs.hpp"
#include "geodesy/geodesy.hpp"
#include "egravity.hpp"
#include <cmath>
#include <array>
#include "gcem.hpp"
#ifdef DEBUG
#  include <limits>
#  include <cfenv>
#  include <cassert>
//#  pragma STDC_FENV_ACCESS on
#endif

namespace {
template <int Nmax> struct IntSquareRoots {
  std::array<double, Nmax + 1> _sqar{};
  constexpr IntSquareRoots() noexcept {
    for (int i = 0; i <= Nmax; i++)
      _sqar[i] = gcem::sqrt(static_cast<double>(i));
  }
  constexpr double operator()(int n) const noexcept {
#ifdef DEBUG
    assert(n >= 0 && n <= Nmax);
#endif
    return _sqar[n];
  }
}; // IntSquareRoots
  inline double Bnm(int n, int m) noexcept {
    const double ar = static_cast<double>( (2*n+1) * (2*n-1) );
    const double pr = static_cast<double>( (n+m)*(n-m) );
    return std::sqrt( ar/pr );
  }
}

int test::gravacc2(const dso::HarmonicCoeffs &cs,
                   const Eigen::Matrix<double, 3, 1> &r, int degree, double Re,
                   double GM, Eigen::Matrix<double, 3, 1> &acc) noexcept {
  const double R = r.norm();
  const double R2 = r.squaredNorm();
  // normalized coordinates
  const double xn = Re * r(0) / R2;
  const double yn = Re * r(1) / R2;
  const double zn = Re * r(2) / R2;
  // factor
  const double rho = Re*Re / R2;

  const int lp_degree = degree + 1; // aka, [0,....degree+1]
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> W(lp_degree+1,
                                                            lp_degree+1);
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> M(lp_degree+1,
                                                            lp_degree+1);
  
  // Look-up table for square roots of integers
  constexpr IntSquareRoots<2*150> isr;
  static_assert(isr(0) == 0e0);

  W.fill_with(0e0);
  M.fill_with(0e0);

  M(0,0) = Re / R;
  W(0,0) = 0e0;
  M(1,0) = isr(3) * zn * M(0,0);
  W(1,0) = 0e0;
  M(1,1) = isr(3) * (xn*M(0,0) - yn*W(0,0));
  W(1,1) = isr(3) * (xn*W(0,0) + yn*M(0,0));

  // first fill column 0;  note that W(n,0) = 0 (already set)
  for (int n=2; n<=lp_degree; n++) {
      const int m = 0;
      M(n, m) = Bnm(n, m) * zn * M(n - 1, m) -
                (Bnm(n, m) / Bnm(n, m - 1)) * rho * M(n - 2, m);
  }
  
  // fill all columns of m >= 1
  for (int m = 1; m <= lp_degree; m++) {
    // M(m,m) and W(m,m)
    if (m>1) {
      M(m,m) = std::sqrt((2*m+1)/m/2e0) * (xn*M(m-1,m-1)-yn*W(m-1,m-1));
      W(m,m) = std::sqrt((2*m+1)/m/2e0) * (xn*W(m-1,m-1)+yn*M(m-1,m-1));
    }

    // warning, if n=m+1 , we do not have a M(n-2,...)
    if (m<lp_degree) {
      const int n = m+1;
      M(n,m) = Bnm(n,m) * zn * M(n-1,m);
      W(n,m) = Bnm(n,m) * zn * W(n-1,m);
    }

    // go on ....
    for (int n = m+2; n <= lp_degree; n++) {
      M(n,m) = Bnm(n,m)*zn*M(n-1,m) - (Bnm(n,m)/Bnm(n,m-1))*rho*M(n-2,m);
      W(n,m) = Bnm(n,m)*zn*W(n-1,m) - (Bnm(n,m)/Bnm(n,m-1))*rho*W(n-2,m);
    }
  }

  // acceleration in cartesian components
  double xacc = 0e0;
  double yacc = 0e0;
  double zacc = 0e0;
  // terms with m = 0
  for (int n = 0; n <= degree; n++) {
    double ar = (2*n+1)*(n+1)*(n+2);
    double pr = 2*(2*(n+1)+1);
    double fac = std::sqrt(ar/pr);
    xacc -= fac * (cs.C(n, 0) * M(n + 1, 1));
    yacc -= fac * (cs.C(n, 0) * W(n + 1, 1));
    zacc -= std::sqrt((2e0*n+1e0)/(2e0*(n+1)+1e0)) * (n + 1) *
            (cs.C(n, 0) * M(n + 1, 0));
  }

  // for m > 0
  for (int m = 1; m <= degree; m++) {
    for (int n = m; n <= degree; n++) {
      double fac1 = ((2e0 * n + 1e0) / (2e0 * (n + 1e0) + 1e0)) * (n + m + 1) *
                    (n + m + 2);
      double fac2 =
          ((m == 1 ? 2e0 : 1e0) * (2e0 * n + 1e0) / (2e0 * (n + 1e0) + 1e0)) *
          (n - m + 2) * (n - m + 1);
      xacc += 0.5 * (std::sqrt(fac1) * (-cs.C(n, m) * M(n + 1, m + 1) +
                                        cs.S(n, m) * W(n + 1, m + 1)) +
                     std::sqrt(fac2) * (-cs.C(n, m) * M(n + 1, m - 1) +
                                        cs.S(n, m) * W(n + 1, m - 1)));

      yacc += .5e0 * (std::sqrt(fac1) * (-cs.C(n, m) * W(n + 1, m + 1) +
                                         cs.S(n, m) * M(n + 1, m + 1)) +
                      std::sqrt(fac2) * (-cs.C(n, m) * W(n + 1, m - 1) +
                                         cs.S(n, m) * W(n + 1, m - 1)));
      fac1 = ((2e0 * n + 1e0) / (2e0 * (n + 1e0) + 1e0)) * (n - m + 1) *
             (n + m + 1);
      zacc += std::sqrt(fac1) *
              (-cs.C(n, m) * M(n + 1, m) - cs.S(n, m) * W(n + 1, m));
    }
  }

  acc = (GM/Re/Re)*Eigen::Matrix<double,3,1>(xacc, yacc,zacc);

  return 0;
}
