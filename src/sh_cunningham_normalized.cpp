#include "cmat2d.hpp"
#include "egravity.hpp"
#include "gcem.hpp"
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

namespace {
// compile-time factors
template <int N> struct NormalizedLegendreFactors {
  std::array<double, N *(N + 1) / 2> fac1{};
  std::array<double, N *(N + 1) / 2> fac2{};
  constexpr int slice(int m) const noexcept { return m * N - m * (m - 1) / 2; }
  constexpr int nm2index(int n, int m) const noexcept {
#ifdef DEBUG
    assert(n > 0 && n < N && m <= n);
#endif
    return slice(m) + (n - m);
  }
  constexpr double &f1(int n, int m) noexcept { return fac1[nm2index(n, m)]; }
  constexpr double f1(int n, int m) const noexcept {
    return fac1[nm2index(n, m)];
  }
  constexpr double &f2(int n, int m) noexcept { return fac2[nm2index(n, m)]; }
  constexpr double f2(int n, int m) const noexcept {
    return fac2[nm2index(n, m)];
  }
  constexpr NormalizedLegendreFactors() noexcept {
    // factors for the recursion ( (2n+1)/2n )^(1/2)
    f1(1, 1) = gcem::sqrt(3e0);
    for (int n = 2; n < N; n++)
      f1(n, n) = gcem::sqrt(static_cast<double>(2 * n + 1) /
                            static_cast<double>(2 * n));

    // factors for the recursion
    for (int m = 0; m < N; m++) {
      for (int n = m + 1; n < N; n++) {
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
} // namespace

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

int test::gravacc2(const dso::HarmonicCoeffs &cs,
                   const Eigen::Matrix<double, 3, 1> &p, int degree, double Re,
                   double GM, Eigen::Matrix<double, 3, 1> &acc) noexcept {

  const int lp_degree = degree + 1; // aka, [0,....degree+1]
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> W(lp_degree + 1,
                                                            lp_degree + 1);
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> M(lp_degree + 1,
                                                            lp_degree + 1);
  W.fill_with(0e0);
  M.fill_with(0e0);

  const NormalizedLegendreFactors<125> F;
  //assert(F.f1(1, 1) == gcem::sqrt(3e0));
  //const Factor F;
  //for (int n=1; n<=lp_degree; n++) {
  //  for (int m=0; m<=n; m++) {
  //    if (F.f1(n,m) != Fct.f1(n,m)) {
  //      printf("ERROR at (%d,%d) runtime gives: %.12f compile-time gives: %.12f diff=%.20e\n", n, m, F.f1(n,m), Fct.f1(n,m), std::abs(F.f1(n,m) - Fct.f1(n,m)));
  //    }
  //  }
  //}

#ifdef KOKO
  {
    const Eigen::Matrix<double, 3, 1> r = p / Re;
    const double rr = std::pow(1 / r.norm(), 2);
    const double x = r.x() * rr;
    const double y = r.y() * rr;
    const double z = r.z() * rr;
    M(0, 0) = 1e280 / r.norm();

    // Recursion diagonal
    // C(n-1,n-1) -> C(n,n)
    for (int n = 1; n <= degree; n++) {
      M(n, n) = F.f1(n, n) * (x * M(n - 1, n - 1) - y * W(n - 1, n - 1));
      W(n, n) = F.f1(n, n) * (y * M(n - 1, n - 1) + x * W(n - 1, n - 1));
    }

    // Recursion secondary diagonal
    // C(n-1,n-1) -> C(n,n-1)
    for (int n = 1; n <= degree; n++) {
      M(n, n - 1) = F.f1(n, n - 1) * z * M(n - 1, n - 1);
      W(n, n - 1) = F.f1(n, n - 1) * z * W(n - 1, n - 1);
    }

    // Recursion others
    // C(n-1,m),C(n-1,m) -> C(n,m)
    if (degree > 1) {
      // for (int m = 0; m < (degree - 1); m++) {
      for (int m = 0; m < lp_degree-1; m++) {
        for (int n = m + 2; n <= degree; n++) {
          M(n, m) =
              F.f1(n, m) * z * M(n - 1, m) + F.f2(n, m) * rr * M(n - 2, m);
          W(n, m) =
              F.f1(n, m) * z * W(n - 1, m) + F.f2(n, m) * rr * W(n - 2, m);
        }
      }
    }

    M.multiply(1e-280);
    W.multiply(1e-280);
  }
#endif

  // Associated Legendre Polynomials
  //dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> W1(lp_degree + 1,
  //                                                          lp_degree + 1);
  //dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> M1(lp_degree + 1,
  //                                                          lp_degree + 1);
  {
    const Eigen::Matrix<double, 3, 1> r = p / Re;
    const double rr = std::pow(1 / r.norm(), 2);
    const double x = r.x() * rr;
    const double y = r.y() * rr;
    const double z = r.z() * rr;

    // start ALF iteration
    M(0, 0) = 1e0 / r.norm();
    M(1, 0) = std::sqrt(3e0)*z*M(0,0);// F.f1(1, 0) * z * M(0, 0);

    // first fill column 0;  note that W(n,0) = 0 (already set)
    for (int n = 2; n <= lp_degree; n++)
      M(n, 0) = F.f1(n, 0) * z * M(n - 1, 0) + F.f2(n, 0) * rr * M(n - 2, 0);

    // fill all columns of m >= 1
    for (int m = 1; m < lp_degree; m++) {
      // M(m,m) and W(m,m) aka, diagonal
      M(m, m) = F.f1(m, m) * (x * M(m - 1, m - 1) - y * W(m - 1, m - 1));
      W(m, m) = F.f1(m, m) * (y * M(m - 1, m - 1) + x * W(m - 1, m - 1));

      // if n=m+1 , we do not have a M(n-2,...) aka sub-diagonal term
      M(m + 1, m) = F.f1(m + 1, m) * z * M(m, m);
      W(m + 1, m) = F.f1(m + 1, m) * z * W(m, m);

      // go on ....
      for (int n = m + 2; n <= lp_degree; n++) {
        M(n, m) = F.f1(n, m) * z * M(n - 1, m) + F.f2(n, m) * rr * M(n - 2, m);
        W(n, m) = F.f1(n, m) * z * W(n - 1, m) + F.f2(n, m) * rr * W(n - 2, m);
      }
    }

    // well, we've left the lst term uncomputed
    {
      const int m = lp_degree;
      M(m, m) = F.f1(m, m) * (x * M(m - 1, m - 1) - y * W(m - 1, m - 1));
      W(m, m) = F.f1(m, m) * (y * M(m - 1, m - 1) + x * W(m - 1, m - 1));
    }

    //M.multiply(1e-280);
    //W.multiply(1e-280);
  } // end computing Legendre
  //int error = 0;
  //for (int m=0; m<=lp_degree; m++) {
  //  for (int n=m+1; n<=lp_degree; n++) {
  //    if (std::abs(M1(n,m) - M(n,m))>1e-16) {
  //      printf("Error M(%d,%d) = %.12e Vs %.12e diff=%.16e\n", n,m,M1(n,m), M(n,m), std::abs(M1(n,m) - M(n,m)));
  //      ++error;
  //    }
  //  }
  //}
  //if (error) return error;

  // acceleration in cartesian components
  acc = Eigen::Matrix<double,3,1>::Zero();
  const int minDegree = 1;
  for(int n=minDegree; n<=degree; n++) {
      // 0. Order
      double wm0 = std::sqrt(static_cast<double>(n+1)*(n+1));
      double wp1 = std::sqrt(static_cast<double>(n+1)*(n+2)) / std::sqrt(2.0);

      double Cm0 = wm0*M(n+1,0);
      double Cp1 = wp1*M(n+1,1); 
      double Sp1 = wp1*W(n+1,1);

      double gx = cs.C(n,0) * (-2*Cp1);
      double gy = cs.C(n,0) * (-2*Sp1);
      double gz = cs.C(n,0) * (-2*Cm0);

      for (int m=1; m<=n; m++) {
        double wm1 = std::sqrt(static_cast<double>(n-m+1)*(n-m+2)) * ((m==1) ? std::sqrt(2.0) : 1.0);
        wm0 = std::sqrt(static_cast<double>(n-m+1)*(n+m+1));
        wp1 = std::sqrt(static_cast<double>(n+m+1)*(n+m+2));

        const double Cm1 = wm1*M(n+1,m-1);
        const double Sm1 = wm1*W(n+1,m-1);
        Cm0 = wm0*M(n+1,m  );
        const double Sm0 = wm0*W(n+1,m  );
        Cp1 = wp1*M(n+1,m+1);
        Sp1 = wp1*W(n+1,m+1);

        gx += cs.C(n,m) * ( Cm1 - Cp1) + cs.S(n,m) * (Sm1 - Sp1);
        gy += cs.C(n,m) * (-Sm1 - Sp1) + cs.S(n,m) * (Cm1 + Cp1);
        gz += cs.C(n,m) * (-2*Cm0)     + cs.S(n,m) * (-2*Sm0);
      }

    acc +=  std::sqrt((2.*n+1.)/(2.*n+3.))*Eigen::Matrix<double,3,1>(gx,gy,gz);
  }

  acc *= GM / (2e0 * Re * Re);


  // terms with m = 0
  //for (int n = minDegree; n <= degree; n++) {
  //  const double wm0 = std::sqrt(static_cast<double>(2 * n + 1) /
  //                               static_cast<double>(2 * n + 3));
  //  const double wp1 =
  //      std::sqrt(static_cast<double>(n + 1) * (n + 2)) / std::sqrt(2e0);
  //  xacc -= wm0 * wp1 * (cs.C(n, 0) * M(n + 1, 1));
  //  yacc -= wm0 * wp1 * (cs.C(n, 0) * W(n + 1, 1));
  //  zacc -= static_cast<double>(n + 1) * (cs.C(n, 0) * M(n + 1, 0));
  //  printf("adding to z(%d,%d) = %.6f * (%.12f * %.12f)\n", n, 0,
  //         static_cast<double>(n + 1), cs.C(n, 0), M(n + 1, 0));
  //}

  //// for m > 0
  //for (int m = 1; m <= degree; m++) {
  //  for (int n = m; n <= degree; n++) {

  //    const double f1 = std::sqrt((n + m + 1) * (n + m + 2) * 1e0);
  //    const double f2 =
  //        std::sqrt((n - m + 1) * (n - m + 2) * (m == 1 ? 2e0 : 1e0));
  //    const double f3 = std::sqrt((n - m + 1) * (n + m + 1) * 1e0);

  //    // TODO need factors dependent on n!
  //    double axt(0e0), ayt(0e0), azt(0e0);

  //    axt =
  //        .5e0 *
  //        (f1 * (-cs.C(n, m) * M(n + 1, m + 1) - cs.S(n, m) * W(n + 1, m + 1)) +
  //         f2 * (cs.C(n, m) * M(n + 1, m - 1) + cs.S(n, m) * W(n + 1, m - 1)));

  //    ayt =
  //        .5e0 *
  //        (f1 * (-cs.C(n, m) * W(n + 1, m + 1) + cs.S(n, m) * M(n + 1, m + 1)) +
  //         f2 * (-cs.C(n, m) * W(n + 1, m - 1) + cs.S(n, m) * W(n + 1, m - 1)));

  //    azt = f3 * (-cs.C(n, m) * M(n + 1, m) - cs.S(n, m) * W(n + 1, m));

  //    const double wm0 = std::sqrt(static_cast<double>(2 * n + 1) /
  //                                 static_cast<double>(2 * n + 3));

  //    xacc += wm0 * axt;
  //    yacc += wm0 * ayt;
  //    zacc += wm0 * azt;
  //  }
  //}

  //acc = (GM / Re / Re) * Eigen::Matrix<double, 3, 1>(xacc, yacc, zacc);

  return 0;
}
