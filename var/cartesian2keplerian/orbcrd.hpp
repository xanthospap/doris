#include "ggeodesy/geoconst.hpp"
#include "ggeodesy/units.hpp"
#include <cmath>

namespace ngpt {
constexpr double GM = 3.9860044e14;
}

namespace math_3d {
inline void cross_product(const double *__restrict__ a,
                          const double *__restrict__ b,
                          double *__restrict__ c) noexcept {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

inline double dot_product(const double *__restrict__ a,
                          const double *__restrict__ b) noexcept {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
}// math_3d

void cartesian2keplerian_tapley(const double *x, const double *xdot) noexcept;
void cartesian2keplerian_schwarz(const double *x, const double *xdot) noexcept;
void cartesian2keplerian_beutler(const double *x, const double *xdot,
                                 double t = 0e0) noexcept;
void cartesian2keplerian_bernese(const double *x, const double *v,
                                 double t = 0e0) noexcept;

/*
inline void rotation_matrix_x(double a, double *R) noexcept {
  const double cosa = std::cos(a);
  const double sina = std::sin(a);
  R[0] = 1e0; R[1] = 0e0;   R[2] = 0e0;
  R[3] = 0e0; R[4] = cosa;  R[5] = sina;
  R[6] = 0e0; R[7] = -sina; R[8] = cosa;
}
inline void rotation_matrix_y(double a, double *R) noexcept {
  const double cosa = std::cos(a);
  const double sina = std::sin(a);
  R[0] = cosa; R[1] = 0e0; R[2] = -sina;
  R[3] = 0e0;  R[4] = 1e0; R[5] = 0e0;
  R[6] = sina; R[7] = 0e0; R[8] = cosa;
}
inline void rotation_matrix_z(double a, double *R) noexcept {
  const double cosa = std::cos(a);
  const double sina = std::sin(a);
  R[0] = cosa; R[1] = sina; R[2] =0e0;
  R[3] = sina; R[4] = cosa; R[5] = 0e0;
  R[6] = 0e0;  R[7] = 0e0;  R[8] = 1e0;
}
inline void matmul3x3(const double *r1, const double* r2, double* r3) noexcept {
  r3[0] = r1[0]*r2[0] + r1[1]*r2[3] + r1[2]*r2[6];
  r3[1] = r1[0]*r2[1] + r1[1]*r2[3] + r1[2]*r2[7];
  r3[2] = r1[0]*r2[2] + r1[1]*r2[5] + r1[2]*r2[8];
  r3[3] = r1[3]*r2[0] + r1[4]*r2[3] + r1[5]*r2[6];
  r3[4] = r1[3]*r2[1] + r1[4]*r2[3] + r1[5]*r2[7];
  r3[5] = r1[3]*r2[2] + r1[4]*r2[5] + r1[5]*r2[8];
  r3[6] = r1[6]*r2[0] + r1[7]*r2[3] + r1[8]*r2[6];
  r3[7] = r1[6]*r2[1] + r1[7]*r2[4] + r1[8]*r2[7];
  r3[8] = r1[6]*r2[2] + r1[7]*r2[5] + r1[8]*r2[8];
  return;
}
inline void matvecmul(const double *m, const double* v, double* res) noexcept {
  res[0] = m[0]*v[0] + m[1]*v[1] + m[2]*v[2];
  res[1] = m[3]*v[0] + m[4]*v[1] + m[5]*v[2];
  res[2] = m[6]*v[0] + m[7]*v[1] + m[8]*v[2];
}
*/
