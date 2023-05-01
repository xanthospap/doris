#include "cmat2d.hpp"
#include "egravity.hpp"
#include "geodesy/geodesy.hpp"
#include "stokes_coeffs.hpp"
#include <array>
#include <cassert>
#include <cmath>
#ifdef DEBUG
#include <cassert>
#include <cfenv>
#include <limits>
#endif

namespace {
struct KahanSum {
  double sum{0e0};
  double c{0e0}; // error
  void operator+(double val) noexcept {
    const double y = val - c;
    const double t = sum + y;
    c = (t - sum) - y;
    sum = t;
  }
  double operator()() const noexcept { return sum; }
}; // KahanSum

/* Max size for ALF factors; if degree is more than this (-2), then it must
 * be augmented. For now, OK
 */
constexpr const int MAX_SIZE_FOR_ALF_FACTORS = 185;

/* Normalized factors for spherical harmonics expansion. These factors
 * depend only on degree/order, not Stokes coeffs or point of expansion
 */
struct NormalizedLegendreFactors {
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> f1;
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> f2;
  std::array<double, MAX_SIZE_FOR_ALF_FACTORS> f3;

  NormalizedLegendreFactors() noexcept
      : f1(MAX_SIZE_FOR_ALF_FACTORS, MAX_SIZE_FOR_ALF_FACTORS),
        f2(MAX_SIZE_FOR_ALF_FACTORS, MAX_SIZE_FOR_ALF_FACTORS) {
    constexpr const int N = MAX_SIZE_FOR_ALF_FACTORS;
    f1.fill_with(0e0);
    f2.fill_with(0e0);

    /* factors for the recursion ( (2n+1)/2n )^(1/2) */
    f1(1, 1) = std::sqrt(3e0);
    for (int n = 2; n < N; n++) {
      f1(n, n) = std::sqrt((2e0 * n + 1e0) / (2e0 * n));
    }

    /* factors for the recursion */
    for (int m = 0; m < N - 1; m++) {
      for (int n = m + 1; n < N; n++) {
        const double f =
            (2e0 * n + 1e0) / static_cast<double>((n + m) * (n - m));
        f1(n, m) = std::sqrt(f * (2e0 * n - 1e0));
        f2(n, m) =
            -std::sqrt(f * (n - m - 1e0) * (n + m - 1e0) / (2e0 * n - 3e0));
      }
    }

    /* factors for acceleration */
    for (int n = 0; n < N; n++)
      f3[n] = std::sqrt((double)(2 * n + 1) / (2 * n + 3));
  }
}; /* NormalizedLegendreFactors */
} /* unnamed namespace */

/*            | dx/dx dx/dy dx/dz |
 * gradient = | dy/dx dy/dy dy/dz |
 *            | dz/dx dz/dy dz/dz |
 */
int gravity_acceleration_impl(
    const dso::StokesCoeffs &cs, const Eigen::Matrix<double, 3, 1> &p,
    int degree, double Re, double GM, Eigen::Matrix<double, 3, 1> &acc,
    Eigen::Matrix<double, 3, 3> &gradient,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &W,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &M) noexcept {

  /* make sure degree is compatible with MAX_SIZE_FOR_ALF_FACTORS */
  if (degree > MAX_SIZE_FOR_ALF_FACTORS - 3) {
    fprintf(stderr,
            "[ERROR] (Static) Size for NormalizedLegendreFactors must be "
            "augmented to perform computation (traceback: %s)\n",
            __func__);
    assert(false);
  }

  /* Factors up to degree/order MAX_SIZE_FOR_ALF_FACTORS; these only depend on
   * degree and order, so they are stored in static storage for subsequent
   * calls
   */
  static const NormalizedLegendreFactors F;

  /* initialize SH coefficients to zero */
  W.fill_with(0e0);
  M.fill_with(0e0);

  /* because we are computing gradient, we need more than degree coefficients.
   * We are actually going for degree+2 coeffs, aka the range [0,degree+2]
   */
  const int lp_degree = degree + 2;

  { /* (kinda) Associated Legendre Polynomials M and W
     * note that since we are going to compute derivatives (of the gravity
     * potential, we will compute M and W up to degree n+2, i.e.
     * M(0,0) to M(lp_degree, lp_degree)
     */
    const Eigen::Matrix<double, 3, 1> r = p / Re;
    const double rr = std::pow(1e0 / r.norm(), 2);
    const double x = r.x() * rr;
    const double y = r.y() * rr;
    const double z = r.z() * rr;

    /* start ALF iteration (can use scaling according to Holmes et al, 2002) */
    constexpr const int __use_factor = true;
    if (__use_factor)
      M(0, 0) = 1e280 / r.norm();
    else
      M(0, 0) = 1e0 / r.norm();

    /* first fill m=0 terms; note that W(n,0) = 0 (already set) */
    M(1, 0) = std::sqrt(3e0) * z * M(0, 0);
    for (int n = 2; n <= lp_degree; n++) {
      M(n, 0) = F.f1(n, 0) * z * M(n - 1, 0) + F.f2(n, 0) * rr * M(n - 2, 0);
    }

    /* fill all elements for order m >= 1 */
    for (int m = 1; m < lp_degree; m++) {
      /* M(m,m) and W(m,m) aka, diagonal */
      M(m, m) = F.f1(m, m) * (x * M(m - 1, m - 1) - y * W(m - 1, m - 1));
      W(m, m) = F.f1(m, m) * (y * M(m - 1, m - 1) + x * W(m - 1, m - 1));

      /* if n=m+1 , we do not have a M(n-2,...) aka sub-diagonal term */
      M(m + 1, m) = F.f1(m + 1, m) * z * M(m, m);
      W(m + 1, m) = F.f1(m + 1, m) * z * W(m, m);

      /* go on .... */
      for (int n = m + 2; n <= lp_degree; n++) {
        M(n, m) = F.f1(n, m) * z * M(n - 1, m) + F.f2(n, m) * rr * M(n - 2, m);
        W(n, m) = F.f1(n, m) * z * W(n - 1, m) + F.f2(n, m) * rr * W(n - 2, m);
      }
    }

    { /* well, we've left the last term uncomputed */
      const int m = lp_degree;
      M(m, m) = F.f1(m, m) * (x * M(m - 1, m - 1) - y * W(m - 1, m - 1));
      W(m, m) = F.f1(m, m) * (y * M(m - 1, m - 1) + x * W(m - 1, m - 1));
    }

    if (__use_factor) {
      M.multiply(1e-280);
      W.multiply(1e-280);
    }
  } /* end computing ALF */

  /* acceleration and gradient in cartesian components */
  acc = Eigen::Matrix<double, 3, 1>::Zero();
  gradient = Eigen::Matrix<double, 3, 3>::Zero();
#ifdef KAHAN_SUM
  KahanSum axs, ays, azs;
#endif

  /* start from smaller terms. note that for degrees m=0,1, we are using
   * seperate loops
   */
  for (int m = degree; m >= 2; --m) {
    for (int n = degree; n >= m; --n) {
      { /* acceleration */
        const double wm1 =
            std::sqrt(static_cast<double>((n - m + 1) * (n - m + 2)));
        const double wm0 =
            std::sqrt(static_cast<double>((n - m + 1) * (n + m + 1)));
        const double wp1 =
            std::sqrt(static_cast<double>((n + m + 1) * (n + m + 2)));

        const double Cm1 = wm1 * M(n + 1, m - 1);
        const double Sm1 = wm1 * W(n + 1, m - 1);
        const double Cm0 = wm0 * M(n + 1, m);
        const double Sm0 = wm0 * W(n + 1, m);
        const double Cp1 = wp1 * M(n + 1, m + 1);
        const double Sp1 = wp1 * W(n + 1, m + 1);

        const double ax = cs.C(n, m) * (Cm1 - Cp1) + cs.S(n, m) * (Sm1 - Sp1);
        const double ay = cs.C(n, m) * (-Sm1 - Sp1) + cs.S(n, m) * (Cm1 + Cp1);
        const double az = cs.C(n, m) * (-2 * Cm0) + cs.S(n, m) * (-2 * Sm0);

#ifdef KAHAN_SUM
        axs + ax *std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
        ays + ay *std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
        azs + az *std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
#endif
        acc += Eigen::Matrix<double, 3, 1>(ax, ay, az) *
               std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
      }
      { /* derivative of acceleration */
        const double wm2 =
            std::sqrt(static_cast<double>((n - m + 1) * (n - m + 2) *
                                          (n - m + 3) * (n - m + 4))) *
            ((m == 2) ? std::sqrt(2.0) : 1.0);
        const double wm1 = std::sqrt(static_cast<double>(
            (n - m + 1) * (n - m + 2) * (n - m + 3) * (n + m + 1)));
        const double wm0 = std::sqrt(static_cast<double>(
            (n - m + 1) * (n - m + 2) * (n + m + 1) * (n + m + 2)));
        const double wp1 = std::sqrt(static_cast<double>(
            (n - m + 1) * (n + m + 1) * (n + m + 2) * (n + m + 3)));
        const double wp2 = std::sqrt(static_cast<double>(
            (n + m + 1) * (n + m + 2) * (n + m + 3) * (n + m + 4)));

        const double Cm2 = wm2 * M(n + 2, m - 2);
        const double Sm2 = wm2 * W(n + 2, m - 2);
        const double Cm1 = wm1 * M(n + 2, m - 1);
        const double Sm1 = wm1 * W(n + 2, m - 1);
        const double Cm0 = wm0 * M(n + 2, m);
        const double Sm0 = wm0 * W(n + 2, m);
        const double Cp1 = wp1 * M(n + 2, m + 1);
        const double Sp1 = wp1 * W(n + 2, m + 1);
        const double Cp2 = wp2 * M(n + 2, m + 2);
        const double Sp2 = wp2 * W(n + 2, m + 2);

        const double gxx = cs.C(n, m) * (Cm2 - 2 * Cm0 + Cp2) +
                           cs.S(n, m) * (Sm2 - 2 * Sm0 + Sp2);
        const double gxy = cs.C(n, m) * (-Sm2 + Sp2) + cs.S(n, m) * (Cm2 - Cp2);
        const double gxz = cs.C(n, m) * (-2 * Cm1 + 2 * Cp1) +
                           cs.S(n, m) * (-2 * Sm1 + 2 * Sp1);
        const double gyy = cs.C(n, m) * (-Cm2 - 2 * Cm0 - Cp2) +
                           cs.S(n, m) * (-Sm2 - 2 * Sm0 - Sp2);
        const double gyz = cs.C(n, m) * (2 * Sm1 + 2 * Sp1) +
                           cs.S(n, m) * (-2 * Cm1 - 2 * Cp1);
        const double gzz = cs.C(n, m) * (4 * Cm0) + cs.S(n, m) * (4 * Sm0);

        gradient += Eigen::Matrix<double, 3, 3>{{gxx, gxy, gxz},
                                                {gxy, gyy, gyz},
                                                {gxz, gyz, gzz}} *
                    std::sqrt((2e0 * n + 1e0) / (2e0 * n + 5e0));
      }
    } /* loop over n */
  }   /* loop over m */

  /* order m = 1
   * for (int n = degree; n >= std::max(1,minDegree); --n)
   * begin summation from smaller terms
   */
  for (int n = degree; n >= 1; --n) {
    const int m = 1;
    { /* acceleration
       * only difference with the generalized formula (aka for random n,m)
       * is in wm1
       */
      const double wm1 =
          std::sqrt(static_cast<double>((n - m + 1) * (n - m + 2))) *
          std::sqrt(2e0);
      const double wm0 =
          std::sqrt(static_cast<double>((n - m + 1) * (n + m + 1)));
      const double wp1 =
          std::sqrt(static_cast<double>((n + m + 1) * (n + m + 2)));

      const double Cm1 = wm1 * M(n + 1, m - 1);
      const double Sm1 = wm1 * W(n + 1, m - 1);
      const double Cm0 = wm0 * M(n + 1, m);
      const double Sm0 = wm0 * W(n + 1, m);
      const double Cp1 = wp1 * M(n + 1, m + 1);
      const double Sp1 = wp1 * W(n + 1, m + 1);

      const double ax = cs.C(n, m) * (Cm1 - Cp1) + cs.S(n, m) * (Sm1 - Sp1);
      const double ay = cs.C(n, m) * (-Sm1 - Sp1) + cs.S(n, m) * (Cm1 + Cp1);
      const double az = cs.C(n, m) * (-2e0 * Cm0) + cs.S(n, m) * (-2 * Sm0);

#ifdef KAHAN_SUM
      axs + ax *std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
      ays + ay *std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
      azs + az *std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
#endif
      acc += Eigen::Matrix<double, 3, 1>(ax, ay, az) *
             std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
    }
    { /* derivative of acceleration */
      const double wm1 =
          std::sqrt(static_cast<double>((n - m + 1) * (n - m + 2) *
                                        (n - m + 3) * (n + m + 1))) *
          std::sqrt(2e0);
      const double wm0 = std::sqrt(static_cast<double>(
          (n - m + 1) * (n - m + 2) * (n + m + 1) * (n + m + 2)));
      const double wp1 = std::sqrt(static_cast<double>(
          (n - m + 1) * (n + m + 1) * (n + m + 2) * (n + m + 3)));
      const double wp2 = std::sqrt(static_cast<double>(
          (n + m + 1) * (n + m + 2) * (n + m + 3) * (n + m + 4)));

      const double Cm1 = wm1 * M(n + 2, m - 1);
      const double Sm1 = wm1 * W(n + 2, m - 1);
      const double Cm0 = wm0 * M(n + 2, m);
      const double Sm0 = wm0 * W(n + 2, m);
      const double Cp1 = wp1 * M(n + 2, m + 1);
      const double Sp1 = wp1 * W(n + 2, m + 1);
      const double Cp2 = wp2 * M(n + 2, m + 2);
      const double Sp2 = wp2 * W(n + 2, m + 2);

      const double gxx =
          cs.C(n, m) * (-3 * Cm0 + Cp2) + cs.S(n, m) * (-Sm0 + Sp2);
      const double gxy = cs.C(n, m) * (-Sm0 + Sp2) + cs.S(n, m) * (-Cm0 - Cp2);
      const double gxz =
          cs.C(n, m) * (-2 * Cm1 + 2 * Cp1) + cs.S(n, m) * (-2 * Sm1 + 2 * Sp1);
      const double gyy =
          cs.C(n, m) * (-Cm0 - Cp2) + cs.S(n, m) * (-3 * Sm0 - Sp2);
      const double gyz =
          cs.C(n, m) * (2 * Sp1) + cs.S(n, m) * (-2 * Cm1 - 2 * Cp1);
      const double gzz = cs.C(n, m) * (4 * Cm0) + cs.S(n, m) * (4 * Sm0);

      gradient += Eigen::Matrix<double, 3, 3>{{gxx, gxy, gxz},
                                              {gxy, gyy, gyz},
                                              {gxz, gyz, gzz}} *
                  std::sqrt((2e0 * n + 1e0) / (2e0 * n + 5e0));
    }
  }

  /* order m = 0
   * begin summation from smaller terms
   */
  for (int n = degree; n >= 0; --n) {
    [[maybe_unused]] const int m = 0;
    { /* acceleration */
      double wm0 = std::sqrt(static_cast<double>(n + 1) * (n + 1));
      double wp1 =
          std::sqrt(static_cast<double>(n + 1) * (n + 2)) / std::sqrt(2e0);

      double Cm0 = wm0 * M(n + 1, 0);
      double Cp1 = wp1 * M(n + 1, 1);
      double Sp1 = wp1 * W(n + 1, 1);

      const double ax = cs.C(n, 0) * (-2e0 * Cp1);
      const double ay = cs.C(n, 0) * (-2e0 * Sp1);
      const double az = cs.C(n, 0) * (-2e0 * Cm0);

#ifdef KAHAN_SUM
      axs + ax *std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
      ays + ay *std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
      azs + az *std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
#endif
      acc += Eigen::Matrix<double, 3, 1>(ax, ay, az) *
             std::sqrt((2e0 * n + 1.) / (2e0 * n + 3e0));
    }
    { /* derivative of acceleration */
      const double wm0 =
          std::sqrt(static_cast<double>((n + 1) * (n + 2) * (n + 1) * (n + 2)));
      const double wp1 = std::sqrt(static_cast<double>((n + 1) * (n + 1) *
                                                       (n + 2) * (n + 3))) /
                         std::sqrt(2e0);
      const double wp2 = std::sqrt(static_cast<double>((n + 1) * (n + 2) *
                                                       (n + 3) * (n + 4))) /
                         std::sqrt(2e0);

      const double Cm0 = wm0 * M(n + 2, 0);
      const double Cp1 = wp1 * M(n + 2, 1);
      const double Sp1 = wp1 * W(n + 2, 1);
      const double Cp2 = wp2 * M(n + 2, 2);
      const double Sp2 = wp2 * W(n + 2, 2);

      const double gxx = cs.C(n, 0) * (-2e0 * Cm0 + 2e0 * Cp2);
      const double gxy = cs.C(n, 0) * (2e0 * Sp2);
      const double gxz = cs.C(n, 0) * (4e0 * Cp1);
      const double gyy = cs.C(n, 0) * (-2e0 * Cm0 - 2e0 * Cp2);
      const double gyz = cs.C(n, 0) * (4e0 * Sp1);
      const double gzz = cs.C(n, 0) * (4e0 * Cm0);

      gradient += Eigen::Matrix<double, 3, 3>{{gxx, gxy, gxz},
                                              {gxy, gyy, gyz},
                                              {gxz, gyz, gzz}} *
                  std::sqrt((2e0 * n + 1e0) / (2e0 * n + 5e0));
    }
  }

  /* scale ... */
  gradient *= GM / (4e0 * Re * Re * Re);
  acc *= GM / (2e0 * Re * Re);
#ifdef KAHAN_SUM
  Eigen::Matrix<double, 3, 1> ks({axs(), ays(), azs()});
  ks *= GM / (2e0 * Re * Re);
  printf("# Acceleration diffs = %.15f %.15f %.15f\n", ks(0) - acc(0),
         ks(1) - acc(1), ks(2) - acc(2));
#endif

  return 0;
}

int dso::gravity_acceleration(const dso::StokesCoeffs &cs,
                              const Eigen::Matrix<double, 3, 1> &p, int degree,
                              double Re, double GM,
                              Eigen::Matrix<double, 3, 1> &acc,
                              Eigen::Matrix<double, 3, 3> &gradient) noexcept {

  const int lp_degree = degree + 3;
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> W(lp_degree,
                                                            lp_degree);
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> M(lp_degree,
                                                            lp_degree);
  return gravity_acceleration_impl(cs, p, degree, Re, GM, acc, gradient, W, M);
}

int dso::gravity_acceleration(
    const dso::StokesCoeffs &cs, const Eigen::Matrix<double, 3, 1> &p,
    int degree, double Re, double GM, Eigen::Matrix<double, 3, 1> &acc,
    Eigen::Matrix<double, 3, 3> &gradient,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *W,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *M) noexcept {

  const int lp_degree = degree + 3;
  if (((W->rows() < lp_degree) || (W->cols() < lp_degree)) ||
      ((M->rows() < lp_degree) || (M->cols() < lp_degree))) {
    fprintf(stderr,
            "[ERROR] Erronous size of M/W coefficient matrices for SH "
            "(traceback: %s)\n",
            __func__);
    return 1;
  }

  return gravity_acceleration_impl(cs, p, degree, Re, GM, acc, gradient, *W,
                                   *M);
}