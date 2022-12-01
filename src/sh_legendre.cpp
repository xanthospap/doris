#include "cmat2d.hpp"
#include "harmonic_coeffs.hpp"
#include "geodesy/geodesy.hpp"
#include <cmath>

namespace {
  inline int kdelta(int i, int j) noexcept { return 1*(i==j); };
  inline double _gnm(int n, int m) noexcept {
    return std::sqrt(static_cast<double>((2 * n + 1) * (2 * n - 1)) /
                     static_cast<int>((n + m) * (n - m)));
  }
  inline double _hnm(int n, int m) noexcept {
    return _gnm(n, m) / _gnm(n - 1, m);
  }
  inline double _knm(int n, int m) noexcept {
    return std::sqrt(
        static_cast<double>((2 - kdelta(0, m)) * (n - m) * (n + m + 1)) / 2e0);
  }
  inline double
  dPdphi(int n, int m, double tanf,
         const dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>
             *const p) noexcept {
    return _knm(n, m) * p->operator()(n, m + 1) -
           m * tanf * p->operator()(n, m);
  }
}

int gravacc1(const dso::HarmonicCoeffs &cs, const Eigen::Matrix<double, 3, 1> &r,
            int degree, double Re) noexcept 
{
  // cartesian to spherical, using the Earth's equatorial radius
  const Eigen::Matrix<double, 3, 1> sph = dso::car2sph(r);

  // geocentric latitude [rad]
  assert(sph(1) != M_PI/2e0);
  const double phi = M_PI/2e0 - sph(1);
  const double cf = std::cos(phi);
  const double sf = std::sin(phi);

  // workspace for Legendre polynomials
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> P(degree + 1,
                                                            degree + 1);

  // compute Legendre polynomials
  P(0,0) = 1e0;
  // compute sectorials, P_nn and subsectorials P(n+1,n)
  for (int n = 1; n < degree; n++) {
    P(n, n) = std::sqrt(static_cast<double>((1 + kdelta(1, n)) * (2 * n + 1)) /
                        (2e0 * n)) *
              cf * P(n - 1, n - 1);
    // subsectorial
    P(n+1,n) = _gnm(n+1,n) * sf * P(n,n) - _hnm(n-1,n) * P(n-2,n);
  }

  // we are now missing the last sectorial
  {
    const int n = degree;
    P(n, n) = std::sqrt(static_cast<double>((1 + kdelta(1, n)) * (2 * n + 1)) /
                        (2e0 * n)) *
              cf * P(n - 1, n - 1);
  }

  // remainding terms, keeping m=const and ascending n
  for (int m=0; m<=degree; m++) {
    for (int n=m+2; n<=degree; n++) {
      P(n,m) = _gnm(n,m) * sf * P(n-1,m) - _hnm(n,m) * P(n-2,m);
    }
  }

  // equation (14) to compute potential V
  const double ar = Re / sph(0);
  const double tanf = std::tan(phi);
  const double sl = std::sin(sph(2));
  const double cl = std::cos(sph(2));
  double Vl = 0e0;
  double Vf = 0e0;
  double Vr = 0e0;
  for (int m=0; m<=degree; m++) {
    double Am0 = 0e0;
    double Bm0 = 0e0;
    double Amf = 0e0;
    double Bmf = 0e0;
    double Amr = 0e0;
    double Bmr = 0e0;
    double arn = 1e0;
    for (int n=m; n<=degree; n++) {
      Am0 += arn * cs.C(n,m) * P(n,m);
      Bm0 += (m==0) ? 0e0 : arn * cs.S(n,m) * P(n,m);
      Amf += arn * cs.C(n,m) * dPdphi(n,m,tanf,&P);
      Bmf += (m==0) ? 0e0 : arn * cs.S(n,m) * dPdphi(n,m,tanf,&P);
      Amr += arn * (n+1) * cs.C(n,m) * P(n,m);
      Bmr += (m==0) ? 0e0 : arn * (n+1) * cs.S(n,m) * P(n,m);
      arn *= ar;
    }
    Vl += m * (Am0 * sl - Bm0 * cl);
    Vf += Amf * cl + Bm0 * sl;
    Vr += Amr * cl + Bmr * sl;
  }

  return 0;
}
