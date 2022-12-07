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
  // compile-time factors
template <int N> struct NormalizedLegendreFactors {
  std::array<double, N *(N + 1) / 2> fac1{};
  std::array<double, N *(N + 1) / 2> fac2{};
  constexpr const int slice(int m) noexcept { return m * N - m * (m - 1) / 2; }
  constexpr const int nm2index(int n, int m) noexcept {
#ifdef DEBUG
    static_assert(n>0 && n<N && m<=n);
#endif
    return slice(m) + (n - m);
  }
  constexpr double &f1(int n, int m) noexcept {
    return fac1[nm2index(n,m)];
  }
  constexpr double f1(int n, int m) const noexcept {
    return fac1[nm2index(n,m)];
  }
  constexpr double &f2(int n, int m) noexcept {
    return fac2[nm2index(n,m)];
  }
  constexpr double f2(int n, int m) const noexcept {
    return fac2[nm2index(n,m)];
  }
  constexpr NormalizedLegendreFactors() noexcept {
    // factors for the recursion ( (2n+1)/2n )^(1/2)
    f1(1, 1) = gcem::sqrt(3e0);
    for (int n = 2; n < N; n++)
      f1(n, n) = gcem::sqrt(static_cast<double>(2 * n + 1) /
                              static_cast<double>(2 * n));

    // factors for the recursion 
    for (int m = 0; m < N; m++) {
      for (int n = m + 1; n < m; n++) {
        double f = static_cast<double>(2 * n + 1) /
                   static_cast<double>((n + m) * (n - m));
        // f1_nm = B_nm
        f1(n, m) = gcem::sqrt(f * (2 * n - 1));
        // f2_nm = B_nm / Bn-1m
        f2(n, m) =
            -gcem::sqrt(f * (n - m - 1e0) * (n + m - 1e0) / (2e0 * n - 3e0));
      }
    }
  }
}; // SHFactors
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
  NormalizedLegendreFactors<150> F;
  static_assert(F.f1(1,1) == gcem::sqrt(3e0));
  
  W.fill_with(0e0);
  M.fill_with(0e0);

  // start ALF iteration
  M(0,0) = Re / R;
  W(0,0) = 0e0;
  M(1,0) = F.f1(0+1,0) * zn * M(0,0);
  
  // first fill column 0;  note that W(n,0) = 0 (already set)
  for (int n=2; n<=lp_degree; n++) {
    M(n,0) = F.f1(n,0) * zn * M(n-1,0) + F.f2(n,0) * rho * M(n-2,0);
    // W(n,0) = F.f1(n,0) * zn * W(n-1,0) + F.f2(n,0) * rho * W(n-2,0);
  }
  
  // fill all columns of m >= 1
  for (int m = 1; m < lp_degree; m++) {
    // M(m,m) and W(m,m) aka, diagonal
    M(m,m) = F.f1(m,m) * (xn * M(m-1,m-1) - yn * W(m-1,m-1));
    W(m,m) = F.f1(m,m) * (yn * M(m-1,m-1) + xn * W(m-1,m-1));

    // if n=m+1 , we do not have a M(n-2,...) aka sub-diagonal term
    M(m+1,m) = F.f1(m+1,m) * zn * M(m,m);
    W(m+1,m) = F.f1(m+1,m) * zn * W(m,m);

    // go on ....
    for (int n = m+2; n <= lp_degree; n++) {
      M(n,m) = F.f1(n,m)*zn*M(n-1,m) - F.f2(n,m)*rho*M(n-2,m);
      W(n,m) = F.f1(n,m)*zn*W(n-1,m) - F.f2(n,m)*rho*W(n-2,m);
    }
  }

  // well, we've left the lst term uncomputed
  {
    const int m = lp_degree;
    M(m,m) = F.f1(m,m) * (xn * M(m-1,m-1) - yn * W(m-1,m-1));
    W(m,m) = F.f1(m,m) * (yn * M(m-1,m-1) + xn * W(m-1,m-1));
  }

  // acceleration in cartesian components
  double xacc = 0e0;
  double yacc = 0e0;
  double zacc = 0e0;
  
  // terms with m = 0
  for (int n = 0; n <= degree; n++) {
    const double wm0 = std::sqrt(static_cast<double>(2 * n + 1) /
                                 static_cast<double>(2 * n + 3));
    const double wp1 =
        std::sqrt(static_cast<double>(n + 1) * (n + 2)) / std::sqrt(2e0);
    xacc -= wm0 * wp1 * (cs.C(n, 0) * M(n + 1, 1));
    yacc -= wm0 * wp1 * (cs.C(n, 0) * W(n + 1, 1));
    zacc -= static_cast<double>(n + 1) * (cs.C(n, 0) * M(n + 1, 0));
  }

  // for m > 0
  for (int m = 1; m <= degree; m++) {
    for (int n = m; n <= degree; n++) {

      //double fac1 = ((2e0 * n + 1e0) / (2e0 * (n + 1e0) + 1e0)) * (n + m + 1) *
      //              (n + m + 2);
      //double fac2 =
      //    ((m == 1 ? 2e0 : 1e0) * (2e0 * n + 1e0) / (2e0 * (n + 1e0) + 1e0)) *
      //    (n - m + 2) * (n - m + 1);
      const double f1 = std::sqrt( (n+m+1)*(n+m+2) * 1e0);
      const double f2 = std::sqrt( (n-m+1)*(n-m+2) * (m==1?2e0:1e0));

      // TODO need factors dependent on n!

      xacc += 0.5e0 * (f1 * (-cs.C(n, m) * M(n + 1, m + 1) -
                                        cs.S(n, m) * W(n + 1, m + 1)) +
                     f2 * (cs.C(n, m) * M(n + 1, m - 1) +
                                        cs.S(n, m) * W(n + 1, m - 1)));

      yacc += .5e0 * (f1 * (-cs.C(n, m) * W(n + 1, m + 1) +
                                         cs.S(n, m) * M(n + 1, m + 1)) +
                      f2 * (-cs.C(n, m) * W(n + 1, m - 1) +
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
