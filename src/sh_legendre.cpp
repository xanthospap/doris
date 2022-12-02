#include "cmat2d.hpp"
#include "harmonic_coeffs.hpp"
#include "geodesy/geodesy.hpp"
#include "egravity.hpp"
#include <cmath>
#ifdef DEBUG
#  include <limits>
#  include <cfenv>
//#  pragma STDC_FENV_ACCESS on
#endif

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

#ifdef DEBUG
int report_fenv(const char *msg) {
  int error = 0;
  int n = std::fetestexcept(FE_ALL_EXCEPT);
  if(n & FE_INVALID) {
    ++error;
    fprintf(stderr, "[ERROR] Floating point error, FE_INVALID raised!\n");
    fprintf(stderr, "%s\n", msg);
  }
  std::feclearexcept(FE_ALL_EXCEPT);
  return error;
}
#endif

int test::gravacc1(const dso::HarmonicCoeffs &cs,
                   const Eigen::Matrix<double, 3, 1> &r, int degree, double Re,
                   Eigen::Matrix<double, 3, 1> &acc) noexcept {
#ifdef DEBUG
  [[maybe_unused]]double snan = std::numeric_limits<double>::signaling_NaN();
  char buf[264];
#endif

  // cartesian to spherical, using the Earth's equatorial radius
  const Eigen::Matrix<double, 3, 1> sph = dso::car2sph(r);

  // geocentric latitude [rad]
  assert(sph(1) != M_PI/2e0);
  const double phi = sph(1);
  const double cf = std::cos(phi);
  const double sf = std::sin(phi);

  // workspace for Legendre polynomials -- note that because we will be
  // computing derivatives, we will need Legendre polynomials for degree
  // + 1
  const int lp_degree = degree + 1; // aka, [0,....degree+1]
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> P(lp_degree+1,
                                                            lp_degree+1);

  // compute Legendre polynomials
  P(0,0) = 1e0;
  P(1,0) = _gnm(1,0) * sf * P(0,0);
  // compute sectorials, P_nn and subsectorials P(n+1,n)
  for (int n = 1; n < lp_degree; n++) {
    P(n, n) = std::sqrt(static_cast<double>((1 + kdelta(1, n)) * (2 * n + 1)) /
                        (2e0 * n)) *
              cf * P(n - 1, n - 1);
#ifdef DEBUG
  sprintf(buf, "on Legendre computation at (%d,%d)", n,n);
  if (report_fenv(buf)) return 1;
#endif
    // subsectorial
    P(n+1,n) = _gnm(n+1,n) * sf * P(n,n); // - _hnm(n-1,n) * P(n-2,n);
#ifdef DEBUG
  sprintf(buf, "on Legendre computation at (%d,%d)", n+1,n);
  if (report_fenv(buf)) return 1;
#endif
  }

  // we are now missing the last sectorial
  {
    const int n = lp_degree;
    P(n, n) = std::sqrt(static_cast<double>((1 + kdelta(1, n)) * (2 * n + 1)) /
                        (2e0 * n)) *
              cf * P(n - 1, n - 1);
  }
#ifdef DEBUG
  if (report_fenv("after Legendre computation (2)")) return 1;
#endif

  // remainding terms, keeping m=const and ascending n
  for (int m=0; m<=lp_degree; m++) {
    for (int n=m+2; n<=lp_degree; n++) {
      P(n,m) = _gnm(n,m) * sf * P(n-1,m) - _hnm(n,m) * P(n-2,m);
    }
  }

#ifdef DEBUG
  if (report_fenv("after Legendre computation (3)")) return 1;
#endif

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
#ifdef DEBUG
  sprintf(buf, "on SH computation at (%d,%d)", n,m);
  if (report_fenv(buf)) return 1;
#endif
    }
    Vl += m * (Am0 * sl - Bm0 * cl);
    Vf += Amf * cl + Bm0 * sl;
    Vr += Amr * cl + Bmr * sl;
  }

  // transform acceleration vector to cartesian coordinates
  // acc = dso::car2sph_rotation_matrix(sph) * Eigen::Matrix<double,3,1>(Vr,Vf,Vl);
  const double _rs  = sph(0);
  // const Eigen::Matrix<double, 3, 1> sph = dso::car2sph(r);
  const Eigen::Matrix<double,3,3> R{
    {cf*cl, -sf*cl/_rs, -sl/(_rs*cf)},
    {cf*sl, -sf*sl/_rs, cl/(_rs*sf)},
    {sf, cf/_rs, 0e0}};
  acc = R * Eigen::Matrix<double,3,1>(Vr,Vf,Vl);

  return 0;
}
