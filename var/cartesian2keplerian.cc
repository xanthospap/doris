#include <cmath>
#include <cassert>
#include <cstdio>
#include "ggeodesy/geoconst.hpp"
#include "ggeodesy/units.hpp"

namespace ngpt {
  constexpr double GM = 3.9860044e14;
}

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

inline void cross_product(const double *a, const double *b, double *c) noexcept {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

inline double norm_squared(const double *v) noexcept {
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

inline double dot_product(const double* a, const double *b) noexcept {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void cartesian2keplerian_2(const double *x, const double *xdot) noexcept {
  // distance
  const double r = std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  // speed
  const double vsq = xdot[0]*xdot[0] + xdot[1]*xdot[1] + xdot[2]*xdot[2];
  const double v = std::sqrt(vsq);
  // radial velocity
  const double ur = (x[0]*xdot[0] + x[1]*xdot[1] + x[2]*xdot[2]) / r;
  // angular momentum
  const double hv[] = {x[1]*xdot[2] - x[2]*xdot[1],
    x[2]*xdot[0] - x[0]*xdot[2], x[0]*xdot[1]-x[1]*xdot[0]};
  // magnitude of angular momentum
  const double h = std::sqrt(hv[0]*hv[0]+hv[1]*hv[1]+hv[2]*hv[2]);
  // inclination, i
  double i = std::acos(hv[2] / h);
  // node line
  const double Nv[] = {-hv[1], hv[0], 0e0};
  // magnitude of node line
  const double N = std::sqrt(Nv[0]*Nv[0] + Nv[1]*Nv[1]);
  // ascending node, Ω
  double Omega = std::acos(Nv[0] / N);
  if (Nv[1] < 0e0) Omega = ngpt::D2PI - Omega;
  // eccentricity vector, e
  double faca = vsq - ngpt::GM / r;
  double facb = r * ur;
  double ev[] = { faca*x[0] - facb*xdot[0], faca*x[1] - facb*xdot[1], faca*x[2] - facb*xdot[2]};
  ev[0] /= ngpt::GM;
  ev[1] /= ngpt::GM;
  ev[2] /= ngpt::GM;
  double e = std::sqrt(ev[0]*ev[0]+ev[1]*ev[1]+ev[2]*ev[2]);
  // argument of perigee, ω
  double omega = std::acos((Nv[0]*ev[0]+Nv[1]*ev[1])/N/e);
  if (ev[2] < 0e0) omega = ngpt::D2PI - omega;
  // true anomaly, θ
  double theta = std::acos((ev[0]*x[0]+ev[1]*x[1]+ev[2]*x[2])/e/r);
  if (ur<0e0) theta = ngpt::D2PI - theta;
  // perigee radii
  // const double rp = ((h*h)/ngpt::GM) * (1e0/(1e0+e));
  // apogee radii
  // const double ra = ((h*h)/ngpt::GM) * (1e0/(1e0-e));
  // semi-major axis, a
  double a = (h*h) / ngpt::GM / (1e0 - e*e);

  // printf("p    = %15.6f\n", p);
  printf("a    = %15.6f\n", a);
  printf("e    = %15.13f\n", e);
  printf("i    = %15.12f\n", ngpt::rad2deg(i));
  printf("Ω    = %15.12f\n", ngpt::rad2deg(Omega));
  printf("ω    = %15.12f\n", ngpt::rad2deg(omega));
  printf("θ    = %15.12f\n", ngpt::rad2deg(theta));
}

void cartesian2keplerian_(const double *x, const double *xdot) noexcept {
  double hv[3];
  cross_product(x, xdot, hv);
  double hsq = dot_product(hv, hv);
  double h = std::sqrt(hsq);

  // ascending node, Ω
  double Omega = std::atan2(hv[0], hv[1]);
  if (Omega<0e0) Omega += ngpt::D2PI;
  
  // inclination, i
  double i = std::acos(hv[2]/h);

  // constant Edot
  double rsq = dot_product(x, x);
  double r = std::sqrt(rsq);
  double Edot = dot_product(xdot, xdot)/2e0 - ngpt::GM / r;

  // argument p
  double p = hsq / ngpt::GM;

  // argument e
  double e = std::sqrt(1e0 + 2e0*hsq*Edot / ngpt::GM / ngpt::GM);

  double Rxi[9], RzO[9], R[9], eov[3];
  rotation_matrix_x(i, Rxi);
  rotation_matrix_z(Omega, RzO);
  matmul3x3(Rxi, RzO, R);
  matvecmul(R, x, eov);
  double u = std::atan2(eov[1], eov[0]);
  if (u<0e0) u += ngpt::D2PI;

  printf("p    = %15.6f\n", p);
  printf("a    = %15.6f\n", 0e0);
  printf("e    = %15.13f\n", e);
  printf("i    = %15.12f\n", ngpt::rad2deg(i));
  printf("Ω    = %15.12f\n", ngpt::rad2deg(Omega));
  printf("ω    = %15.12f\n", ngpt::rad2deg(0e0));
}

/// https://downloads.rene-schwarz.com/download/M002-Cartesian_State_Vectors_to_Keplerian_Orbit_Elements.pdf
void cartesian2keplerian(const double *x, const double *v) noexcept {
  
  // Calculate orbital momentum vector
  double hv[3];
  cross_product(x, v, hv);
  const double hsq = dot_product(hv,hv);
  const double h = std::sqrt(hsq);
  const double rsq = dot_product(x, x);
  const double r = std::sqrt(rsq);

  // eccentricity vector
  double ev[3];
  cross_product(v, hv, ev);
  ev[0] /= ngpt::GM;
  ev[1] /= ngpt::GM;
  ev[2] /= ngpt::GM;
  ev[0] -= (x[0] / r);
  ev[1] -= (x[1] / r);
  ev[2] -= (x[2] / r);

  // Determine the vector n pointing towards the ascending node
  double n[] = { -hv[1], hv[0], 0e0};

  //  true anomaly
  double ta = std::acos(dot_product(ev, x)/ e / r);

  // Orbital inclination, i
  double i = std::acos(h[2] / h_norm);
  assert(i>=0e0 && i<= ngpt::DPI); // not needed; acos is guranteed (if no 
  // errors occur) to return a value in [0,π]

  // longtitude or right ascension of the acending node, Ω
  const double h_xy = std::sqrt(h[0]*h[0] + h[1]*h[1]);
  double sinOmega = h[0] / h_xy;
  double cosOmega = -h[1] / h_xy;
  double Omega = std::atan2(sinOmega, cosOmega);
  if (Omega<0e0) Omega += ngpt::D2PI;
  assert(Omega>=0e0 && Omega<=ngpt::D2PI);

  // (intermidiate quantity) r0
  double r0 = std::sqrt(norm_squared(x));

  // energy per unit mass
  double ksi = norm_squared(v) / 2e0 - ngpt::GM / r0;

  // semi-major axis, a
  double a = -ngpt::GM / (2e0*ksi);

  // eccentricity of the orbit, e
  double e = std::sqrt(1e0 + (2e0*ksi*h_sq)/ngpt::GM/ngpt::GM);

  // semilatus rectum, p
  double p = h_sq / ngpt::GM;

  // true anomaly, ω
  double cosf = (p-r0)/(e*r0);
  double sinf = (p/h_norm*e)*((x[0]*v[0]+x[1]*v[1]+x[2]*v[2])/r0);
  double cosof = (x[0]/r0)*std::cos(Omega) + (x[1]/r0)*std::sin(Omega);
  double sinof = x[2] / (r0*std::sin(i));
  double omega_opf = std::atan2(sinof, cosof);
  if (omega_opf<0e0) omega_opf += ngpt::D2PI;
  double omega_f = std::atan2(sinf, cosf);
  if (omega_f<0e0) omega_f += ngpt::D2PI;
  double omega = omega_opf - omega_f;
  if (omega<0e0) omega += ngpt::D2PI;

  /*
  // eccentric anomaly, E0
  double cosE0 = (r0/a)*cosf + e;
  double sinE0 = (r0/b)*sinf;
  double E0 = std::atan2(sinE0, cosE0);
  if (E0<0e0) E0 += ngpt::D2PI;

  // mean anomaly, M0
  double M0 = E0 - e*std::sin(E0);
  */

  printf("p    = %15.6f\n", p);
  printf("a    = %15.6f\n", a);
  printf("e    = %15.13f\n", e);
  printf("i    = %15.12f\n", ngpt::rad2deg(i));
  printf("Ω    = %15.12f\n", ngpt::rad2deg(Omega));
  printf("ω    = %15.12f\n", ngpt::rad2deg(omega));
  //printf("M0   = %15.12f\n", M0);
} 

int main() {
  double x[] = {5492000.34e0,
                      3984001.40e0,
                         2955.81e0};
  double v[] = {  -3931.046491e0,
                         5498.676921,
                         3665.980697};

  printf("/*-----------------------------------------------------------------*/\n");
  cartesian2keplerian(x, v);
  printf("/*-----------------------------------------------------------------*/\n");
  cartesian2keplerian_(x, v);
  printf("/*-----------------------------------------------------------------*/\n");
  cartesian2keplerian_2(x, v);
  
  double x2[] = {-6045000.0e0,
                -3490000.0e0,
                +2500000.0e0};
  double v2[] = {-3457.0e0,
                +6618.0e0,
                 2533.0e0};
  
  printf("/*-----------------------------------------------------------------*/\n");
  cartesian2keplerian(x2, v2);
  printf("/*-----------------------------------------------------------------*/\n");
  cartesian2keplerian_(x2, v2);
  printf("/*-----------------------------------------------------------------*/\n");
  cartesian2keplerian_2(x2, v2);
  return 0;
}
