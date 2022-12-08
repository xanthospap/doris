#include "cmat2d.hpp"
#include "egravity.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "harmonic_coeffs.hpp"
#include <cmath>
#ifdef DEBUG
#include <cfenv>
#include <limits>
#include "geodesy/units.hpp"
//#  pragma STDC_FENV_ACCESS on
#endif

namespace {
inline int kdelta(int i, int j) noexcept { return (i == j); };
inline double _gnm(int n, int m) noexcept {
#ifdef DEBUG
  assert(n != m);
#endif
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
  return _knm(n, m) * p->operator()(n, m + 1) - m * tanf * p->operator()(n, m);
}
inline double A(int n, int m) noexcept {
  const double u = (2 * n - 1) * (2 * n + 1);
  const double d = (n - m) * (n + m);
  return std::sqrt(u / d);
}
inline double B(int n, int m) noexcept {
  const double u = (2 * n + 1) * (n + m - 1) * (n - m - 1);
  const double d = (n - m) * (n + m) * (2 * n - 3);
  return std::sqrt(u / d);
}
} // namespace

#ifdef DEBUG
int report_fenv(const char *msg) {
  int error = 0;
  int n = std::fetestexcept(FE_ALL_EXCEPT);
  if (n & FE_INVALID) {
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
                   double GM, Eigen::Matrix<double, 3, 1> &acc) noexcept {
#ifdef DEBUG
  [[maybe_unused]] double snan = std::numeric_limits<double>::signaling_NaN();
  char buf[264];
#endif

  // cartesian to spherical, using the Earth's equatorial radius
  const Eigen::Matrix<double, 3, 1> sph = dso::car2sph(r);
  printf("Spherical coordinates are: %.12f %.12f %.12f\n", sph(0), dso::rad2deg(sph(1)), dso::rad2deg(sph(2)));

  // geocentric latitude [rad]
  assert(sph(1) != M_PI / 2e0);
  //
  //double theta = std::atan(r(0)*r(0)+r(1)*r(1));
  //if (r(2)<0e0) theta = M_PI + theta;
  //const double phi = M_PI / 2e0 - theta;
  //
  const double phi = sph(1);
  const double cf = std::cos(phi);
  const double sf = std::sin(phi);

  // workspace for Legendre polynomials -- note that because we will be
  // computing derivatives, we will need Legendre polynomials for degree
  // + 1
  const int lp_degree = degree + 1; // aka, [0,....degree+1]
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> P(lp_degree + 1,
                                                            lp_degree + 1);

  // Associated Legendre Polynomials, normalized
  P(0, 0) = 1e0;
  P(1, 0) = _gnm(1, 0) * sf * P(0, 0);
  for (int m = 0; m <= lp_degree; m++) {
    // handling column m (term under sectorial)
    if (m > 0 && m < lp_degree)
      P(m + 1, m) = _gnm(m + 1, m) * sf * P(m, m);
    // all terms (of same column m)
    for (int n = m + 2; n <= lp_degree; n++)
      P(n, m) = _gnm(n, m) * sf * P(n - 1, m) - _hnm(n, m) * P(n - 2, m);
    // sectorial term for next column
    if (m < lp_degree) {
      const int k = m + 1;
      P(k, k) =
          std::sqrt((1e0 + (k == 1)) * (2 * k + 1) / 2e0 / k) * cf * P(m, m);
    }
  }

#ifdef DEBUG
  if (report_fenv("after Legendre computation (3)"))
    return 1;
#endif

  // equation (14) to compute potential V
  const double ar = Re / sph(0);
  const double tanf = std::tan(phi);
  const double sl = std::sin(sph(2));
  const double cl = std::cos(sph(2));
  double Vl = 0e0;
  double Vf = 0e0;
  double Vr = 0e0;
  for (int m = 0; m <= degree; m++) {
    double Am0 = 0e0;
    double Bm0 = 0e0;
    double Amf = 0e0;
    double Bmf = 0e0;
    double Amr = 0e0;
    double Bmr = 0e0;
    double arn = (m == 0) ? 1e0 : std::pow(ar, m);
    for (int n = m; n <= degree; n++) {
      Am0 += arn * cs.C(n, m) * P(n, m);
      Bm0 += (m == 0) ? 0e0 : arn * cs.S(n, m) * P(n, m);
      Amf += arn * cs.C(n, m) * dPdphi(n, m, tanf, &P);
      Bmf += (m == 0) ? 0e0 : arn * cs.S(n, m) * dPdphi(n, m, tanf, &P);
      Amr += arn * (n + 1) * cs.C(n, m) * P(n, m);
      Bmr += (m == 0) ? 0e0 : arn * (n + 1) * cs.S(n, m) * P(n, m);
      arn *= ar;
#ifdef DEBUG
      sprintf(buf, "on SH computation at (%d,%d)", n, m);
      if (report_fenv(buf))
        return 1;
#endif
    } // end loop on m
    Vl += m * (Am0 * sl - Bm0 * cl);
    Vf += Amf * cl + Bm0 * sl;
    Vr += Amr * cl + Bmr * sl;
  } // end loop on n
  Vl *= -GM / r.norm();
  Vf *= GM / r.norm();
  Vr *= -GM / r.squaredNorm();

  // transform acceleration vector to cartesian coordinates
  const double x = r(0);
  const double y = r(1);
  const double z = r(2);
  const double xy = std::sqrt(x * x + y * y);
  const double xy2 = x * x + y * y;
  const double R = sph(0);
  const double R2 = R * R;
  acc = Eigen::Matrix<double, 3, 1>(
      (Vr / R - (z / xy) * (Vf / R2)) * x - (y / xy2) * Vl,
      (Vr / R - (z / xy) * (Vf / R2)) * y + (x / xy2) * Vl,
      (Vr / R) * z + (xy / R2) * Vf);

  return 0;
}
