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
#include <functional>
#include <thread>
#include <mutex>

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
    assert(n >= 0 && n < N);
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
        // double f = static_cast<double>(2 * n + 1) /
        //           static_cast<double>((n + m) * (n - m));
        //// f1_nm = B_nm
        // f1(n, m) = std::sqrt(f * (2 * n - 1));
        //// f2_nm = B_nm / Bn-1m
        // f2(n, m) =
        //     - std::sqrt(f * (n - m - 1e0) * (n + m - 1e0) / (2e0 * n - 3e0));
        //  -- alternative ... keep integers as late as possible
        const long n1 = 2 * n + 1;
        const long n2 = (n + m) * (n - m);
        const long f11 = n1 * (2 * n - 1);
        const long f12 = n2;
        f1(n, m) = std::sqrt((double)f11 / (double)f12);

        const long f21 = n1 * (n - m - 1) * (n + m - 1);
        const long f22 = n2 * (2 * n - 3);
        f2(n, m) = -std::sqrt((double)f21 / (double)f22);
      }
    }

    // factors for acceleration
    for (int n = 0; n < N; n++)
      fac3[n] = std::sqrt((double)(2 * n + 1) / (double)(2 * n + 3));
  }
}; // SHFactors
} // namespace
void order0(int degree, const dso::HarmonicCoeffs &cs,
           const dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &M,
           const dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &W,
           Eigen::Matrix<double, 3, 1> &acc,
           Eigen::Matrix<double, 3, 3> &gradient, int minDegree) noexcept;
void order1(int degree, const dso::HarmonicCoeffs &cs,
           const dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &M,
           const dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &W,
           Eigen::Matrix<double, 3, 1> &acc,
           Eigen::Matrix<double, 3, 3> &gradient) noexcept;
void ordern(int degree, int minm, int maxm, const dso::HarmonicCoeffs &cs,
           const dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &M,
           const dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &W,
           Eigen::Matrix<double, 3, 1> &acc,
           Eigen::Matrix<double, 3, 3> &gradient) noexcept;

static const NormalizedLegendreFactors<125> F;
std::mutex AccMtx;
std::mutex GrdMtx;

/*            | dx/dx dx/dy dx/dz |
 * gradient = | dy/dx dy/dy dy/dz |
 *            | dz/dx dz/dy dz/dz |
 */
int test::gravacc_prl(const dso::HarmonicCoeffs &cs,
                   const Eigen::Matrix<double, 3, 1> &p, int degree, double Re,
                   double GM, Eigen::Matrix<double, 3, 1> &acc,
                   Eigen::Matrix<double, 3, 3> &gradient) noexcept {

  const int lp_degree = degree + 2; // aka, [0,....degree+1]
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> W(lp_degree + 1,
                                                            lp_degree + 1);
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> M(lp_degree + 1,
                                                            lp_degree + 1);
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
    M(1, 0) = std::sqrt(3e0) * z * M(0, 0); // F.f1(1, 0) * z * M(0, 0);

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

  // acceleration and gradient in cartesian components
  acc = Eigen::Matrix<double, 3, 1>::Zero();
  gradient = Eigen::Matrix<double,3,3>::Zero();
  [[maybe_unused]] const int minDegree = 1;



  {
    int num_threads = 4;
    int tcount = (degree - 2) / num_threads;
    std::vector<std::thread> tvec; tvec.reserve(6);
    //for (int i=0; i<tcount; i++) {
      tvec.emplace_back(ordern, degree, 2, 2 + tcount - 1, std::cref(cs), std::cref(M),
                        std::cref(W), std::ref(acc), std::ref(gradient));
      tvec.emplace_back(ordern, degree, 2 + tcount, 2 + 2 * tcount - 1, std::cref(cs),
                        std::cref(M), std::cref(W), std::ref(acc),
                        std::ref(gradient));
      tvec.emplace_back(ordern, degree, 2 + 2 * tcount, 2 + 3 * tcount - 1,
                        std::cref(cs), std::cref(M), std::cref(W),
                        std::ref(acc), std::ref(gradient));
      tvec.emplace_back(ordern, degree, 2 + 3 * tcount, degree, std::cref(cs),
                        std::cref(M), std::cref(W), std::ref(acc),
                        std::ref(gradient));
      tvec.emplace_back(order1, degree, std::cref(cs), std::cref(M),
                        std::cref(W), std::ref(acc), std::ref(gradient));
      tvec.emplace_back(order0, degree, std::cref(cs), std::cref(M),
                        std::cref(W), std::ref(acc), std::ref(gradient),
                        minDegree);
      //std::thread t0(order0, degree, std::cref(cs), std::cref(M), std::cref(W),
      //               std::ref(acc), std::ref(gradient), minDegree);
      //std::thread t1(order1, degree, std::cref(cs), std::cref(M), std::cref(W),
      //               std::ref(acc), std::ref(gradient));
      //std::thread tm(ordern, degree, 2, degree, std::cref(cs), std::cref(M),
      //               std::cref(W), std::ref(acc), std::ref(gradient));
    //}
    for (auto &t : tvec) t.join();
  }

  // scale ...
  gradient *= GM/(4e0*Re*Re*Re);
  acc *= GM/(2e0*Re*Re);

  return 0;
}

void order0(int degree, const dso::HarmonicCoeffs &cs,
           const dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &M,
           const dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &W,
           Eigen::Matrix<double, 3, 1> &acc,
           Eigen::Matrix<double, 3, 3> &gradient, int minDegree) noexcept {
  // order m = 0
  for (int n = degree; n >= minDegree;
       --n) { // begin summation from smaller terms
    [[maybe_unused]] const int m = 0;

    // acceleration
    {
      double wm0 = std::sqrt(static_cast<double>(n + 1) * (n + 1));
      double wp1 =
          std::sqrt(static_cast<double>(n + 1) * (n + 2)) / std::sqrt(2e0);

      double Cm0 = wm0 * M(n + 1, 0);
      double Cp1 = wp1 * M(n + 1, 1);
      double Sp1 = wp1 * W(n + 1, 1);

      const double ax = cs.C(n, 0) * (-2e0 * Cp1);
      const double ay = cs.C(n, 0) * (-2e0 * Sp1);
      const double az = cs.C(n, 0) * (-2e0 * Cm0);

      std::scoped_lock lock(AccMtx);     
      acc += Eigen::Matrix<double, 3, 1>(ax, ay, az) *
             std::sqrt((2e0 * n + 1.) / (2e0 * n + 3e0));
      //const auto tmp = Eigen::Matrix<double, 3, 1>(ax, ay, az);
      //printf("Contribution of order %3d is %.12f\n", 0, tmp.norm());
    }

    // derivative of acceleration
    {
      const double wm0 =
          std::sqrt(static_cast<double>(n + 1) * (n + 2) * (n + 1) * (n + 2));
      const double wp1 =
          std::sqrt(static_cast<double>(n + 1) * (n + 1) * (n + 2) * (n + 3)) /
          std::sqrt(2e0);
      const double wp2 =
          std::sqrt(static_cast<double>(n + 1) * (n + 2) * (n + 3) * (n + 4)) /
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

      std::scoped_lock lock(GrdMtx);     
      gradient += Eigen::Matrix<double, 3, 3>{{gxx, gxy, gxz},
                                             {gxy, gyy, gyz},
                                             {gxz, gyz, gzz}} *
                 std::sqrt((2e0 * n + 1e0) / (2e0 * n + 5e0));
    }
  }
  return;
}

void order1(int degree, const dso::HarmonicCoeffs &cs,
           const dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &M,
           const dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &W,
           Eigen::Matrix<double, 3, 1> &acc,
           Eigen::Matrix<double, 3, 3> &gradient) noexcept {
  // order m = 1
  for (int n = degree; n >= 1; --n) { // begin summation from smaller terms
    const int m = 1;
    // acceleration
    {
      // only difference with the generalized formula (aka for random n,m)
      // is in wm1
      const double wm1 =
          std::sqrt(static_cast<double>(n - m + 1) * (n - m + 2)) *
          std::sqrt(2e0);
      const double wm0 =
          std::sqrt(static_cast<double>(n - m + 1) * (n + m + 1));
      const double wp1 =
          std::sqrt(static_cast<double>(n + m + 1) * (n + m + 2));

      const double Cm1 = wm1 * M(n + 1, m - 1);
      const double Sm1 = wm1 * W(n + 1, m - 1);
      const double Cm0 = wm0 * M(n + 1, m);
      const double Sm0 = wm0 * W(n + 1, m);
      const double Cp1 = wp1 * M(n + 1, m + 1);
      const double Sp1 = wp1 * W(n + 1, m + 1);

      const double ax = cs.C(n, m) * (Cm1 - Cp1) + cs.S(n, m) * (Sm1 - Sp1);
      const double ay = cs.C(n, m) * (-Sm1 - Sp1) + cs.S(n, m) * (Cm1 + Cp1);
      const double az = cs.C(n, m) * (-2e0 * Cm0) + cs.S(n, m) * (-2 * Sm0);

      std::scoped_lock lock(AccMtx);     
      acc += Eigen::Matrix<double, 3, 1>(ax, ay, az) *
            std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
      //const auto tmp = Eigen::Matrix<double, 3, 1>(ax, ay, az);
      //printf("Contribution of order %3d is %.12f\n", 1, tmp.norm());
    }
    // derivative of acceleration
    {
      const double wm1 = std::sqrt(static_cast<double>(n - m + 1) *
                                   (n - m + 2) * (n - m + 3) * (n + m + 1)) *
                         std::sqrt(2e0);
      const double wm0 = std::sqrt(static_cast<double>(n - m + 1) *
                                   (n - m + 2) * (n + m + 1) * (n + m + 2));
      const double wp1 = std::sqrt(static_cast<double>(n - m + 1) *
                                   (n + m + 1) * (n + m + 2) * (n + m + 3));
      const double wp2 = std::sqrt(static_cast<double>(n + m + 1) *
                                   (n + m + 2) * (n + m + 3) * (n + m + 4));

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

      std::scoped_lock lock(GrdMtx);     
      gradient += Eigen::Matrix<double, 3, 3>{{gxx, gxy, gxz},
                                             {gxy, gyy, gyz},
                                             {gxz, gyz, gzz}} *
                 std::sqrt((2e0 * n + 1e0) / (2e0 * n + 5e0));
    }
  }
  return;
}

void ordern(int degree, int minm, int maxm, const dso::HarmonicCoeffs &cs,
           const dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &M,
           const dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &W,
           Eigen::Matrix<double, 3, 1> &acc,
           Eigen::Matrix<double, 3, 3> &gradient) noexcept {
  // start from smaller terms. note that for degrees m=0,1, we are using
  // seperate loops
  for (int m = maxm; m >= std::min(2,minm); --m) {
  for (int n = degree; n >= m; --n) {
    // acceleration
    {
      const double wm1 =
          std::sqrt(static_cast<double>(n - m + 1) * (n - m + 2));
      const double wm0 =
          std::sqrt(static_cast<double>(n - m + 1) * (n + m + 1));
      const double wp1 =
          std::sqrt(static_cast<double>(n + m + 1) * (n + m + 2));

      const double Cm1 = wm1 * M(n + 1, m - 1);
      const double Sm1 = wm1 * W(n + 1, m - 1);
      const double Cm0 = wm0 * M(n + 1, m);
      const double Sm0 = wm0 * W(n + 1, m);
      const double Cp1 = wp1 * M(n + 1, m + 1);
      const double Sp1 = wp1 * W(n + 1, m + 1);

      const double ax = cs.C(n, m) * (Cm1 - Cp1) + cs.S(n, m) * (Sm1 - Sp1);
      const double ay = cs.C(n, m) * (-Sm1 - Sp1) + cs.S(n, m) * (Cm1 + Cp1);
      const double az = cs.C(n, m) * (-2 * Cm0) + cs.S(n, m) * (-2 * Sm0);

      std::scoped_lock lock(AccMtx);     
      acc += Eigen::Matrix<double, 3, 1>(ax, ay, az) *
             std::sqrt((2e0 * n + 1e0) / (2e0 * n + 3e0));
      // const auto tmp = Eigen::Matrix<double, 3, 1>(ax, ay, az);
      // printf("Contribution of order %3d is %.12f\n", m, tmp.norm());

    }
    // derivative of acceleration
    {
      const double wm2 = std::sqrt(static_cast<double>(n - m + 1) *
                                   (n - m + 2) * (n - m + 3) * (n - m + 4)) *
                         ((m == 2) ? std::sqrt(2.0) : 1.0);
      const double wm1 = std::sqrt(static_cast<double>(n - m + 1) *
                                   (n - m + 2) * (n - m + 3) * (n + m + 1));
      const double wm0 = std::sqrt(static_cast<double>(n - m + 1) *
                                   (n - m + 2) * (n + m + 1) * (n + m + 2));
      const double wp1 = std::sqrt(static_cast<double>(n - m + 1) *
                                   (n + m + 1) * (n + m + 2) * (n + m + 3));
      const double wp2 = std::sqrt(static_cast<double>(n + m + 1) *
                                   (n + m + 2) * (n + m + 3) * (n + m + 4));

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
      const double gxz =
          cs.C(n, m) * (-2 * Cm1 + 2 * Cp1) + cs.S(n, m) * (-2 * Sm1 + 2 * Sp1);
      const double gyy = cs.C(n, m) * (-Cm2 - 2 * Cm0 - Cp2) +
                         cs.S(n, m) * (-Sm2 - 2 * Sm0 - Sp2);
      const double gyz =
          cs.C(n, m) * (2 * Sm1 + 2 * Sp1) + cs.S(n, m) * (-2 * Cm1 - 2 * Cp1);
      const double gzz = cs.C(n, m) * (4 * Cm0) + cs.S(n, m) * (4 * Sm0);

      std::scoped_lock lock(GrdMtx);     
      gradient += Eigen::Matrix<double, 3, 3>{{gxx, gxy, gxz},
                                              {gxy, gyy, gyz},
                                              {gxz, gyz, gzz}} *
                  std::sqrt((2e0 * n + 1e0) / (2e0 * n + 5e0));
    }
  } // loop over n
  } // loop over m

  return;
}
