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

namespace {
// compile-time factors
template <int N> struct NormalizedLegendreFactors {
  std::array<double, N *(N + 1) / 2> fac1{};
  std::array<double, N *(N + 1) / 2> fac2{};
  std::array<double, N> fac3{};
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
  constexpr double f3(int n) const noexcept {
#ifdef DEBUG
    assert( n>=0 && n< N );
#endif
     return fac3[n];
  }
  constexpr NormalizedLegendreFactors() noexcept {
    // factors for the recursion ( (2n+1)/2n )^(1/2)
    f1(1, 1) = std::sqrt(3e0);
    for (int n = 2; n < N; n++)
      f1(n, n) = std::sqrt(static_cast<double>(2 * n + 1) /
                            static_cast<double>(2 * n));

    // factors for the recursion
    for (int m = 0; m < N; m++) {
      for (int n = m + 1; n < N; n++) {
        // --- a little easier .....
        //double f = static_cast<double>(2 * n + 1) /
        //           static_cast<double>((n + m) * (n - m));
        //// f1_nm = B_nm
        //f1(n, m) = std::sqrt(f * (2 * n - 1));
        //// f2_nm = B_nm / Bn-1m
        //f2(n, m) =
        //    - std::sqrt(f * (n - m - 1e0) * (n + m - 1e0) / (2e0 * n - 3e0));
        // -- alternative ... keep integers as late as possible
        const long n1 = 2 * n + 1;
        const long n2 = (n + m) * (n - m);
        const long f11 = n1 * (2*n-1);
        const long f12 = n2;
        f1(n,m) = std::sqrt( (double)f11 / (double)f12 );

        const long f21 = n1 * (n - m - 1) * (n + m - 1);
        const long f22 = n2 * (2 * n - 3);
        f2(n,m) = -std::sqrt( (double)f21 / (double)f22 );
      }
    }

    // factors for acceleration
    for (int n=0; n<N; n++)
      fac3[n] = std::sqrt( (double)(2*n+1) / (double)(2*n+3) );
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

static const NormalizedLegendreFactors<125> F;

int test::gravacc2(const dso::HarmonicCoeffs &cs,
                   const Eigen::Matrix<double, 3, 1> &p, int degree, double Re,
                   double GM, Eigen::Matrix<double, 3, 1> &acc) noexcept {

  const int lp_degree = degree + 1; // aka, [0,....degree+1]
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> W(lp_degree + 2,
                                                            lp_degree + 2);
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> M(lp_degree + 2,
                                                            lp_degree + 2);
  W.fill_with(0e0);
  M.fill_with(0e0);

  // Associated Legendre Polynomials
  {
    const Eigen::Matrix<double, 3, 1> r = p / Re;
    const double rr = std::pow(1 / r.norm(), 2);
    const double x = r.x() * rr;
    const double y = r.y() * rr;
    const double z = r.z() * rr;

    // start ALF iteration
    const int __use_factor = false;
    if (__use_factor)
    M(0, 0) = 1e280 / r.norm();
    else
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

    if (__use_factor) {
    M.multiply(1e-280);
    W.multiply(1e-280);
    }
  } // end computing Legendre

  // acceleration in cartesian components
  acc = Eigen::Matrix<double,3,1>::Zero();
  [[maybe_unused]]const int minDegree = 2;

  for (int m = degree; m >= 1; --m) {
    double acx = 0e0, acy = 0e0, acz = 0e0;
    for (int n = degree; n >= m; --n) {
      const long arxy1 = (n + m + 1) * (n + m + 2);
      const long arxy2 = (n - m + 1) * (n - m + 2) * ((m == 1)?2:1);
      const long arzz1 = (n + m + 1) * (n - m + 1);

      acx +=
          F.f3(n) *
          (std::sqrt((double)arxy1) *
               (-cs.C(n, m) * M(n + 1, m + 1) - cs.S(n, m) * W(n + 1, m + 1)) +
           std::sqrt((double)arxy2) *
               (cs.C(n, m) * M(n + 1, m - 1) + cs.S(n, m) * W(n + 1, m - 1)));
      acy +=
          F.f3(n) *
          (std::sqrt((double)arxy1) *
               (-cs.C(n, m) * W(n + 1, m + 1) + cs.S(n, m) * M(n + 1, m + 1)) +
           std::sqrt((double)arxy2) *
               (-cs.C(n, m) * W(n + 1, m - 1) + cs.S(n, m) * M(n + 1, m - 1)));
      
      acz += 2e0*F.f3(n) * std::sqrt((double)arzz1) *
             (-cs.C(n, m) * M(n + 1, m) - cs.S(n, m) * W(n + 1, m));
    
    } // loop over n
    acc += Eigen::Matrix<double, 3, 1>(acx, acy, acz);
  } // loop over m
 
  // order m = 0
  double gxt = 0e0, gyt = 0e0, gzt = 0e0;
  for (int n=degree; n>=minDegree; --n) { // begin summation from smaller terms
    // const int m = 0;
    const long n3 = (n+1)*(n+2);

    const double fgx = (double)n3 / 2e0;
    gxt += cs.C(n, 0) * F.f3(n) * M(n + 1, 1) * std::sqrt(fgx);
    gyt += cs.C(n, 0) * F.f3(n) * W(n + 1, 1) * std::sqrt(fgx);
    gzt += cs.C(n, 0) * F.f3(n) * M(n + 1, 0) * (double)(n + 1);
  }
  
  acc += (-2e0 * Eigen::Matrix<double, 3, 1>(gxt, gyt, gzt));
  acc *= GM / (2 * Re * Re);

  return 0;
}
