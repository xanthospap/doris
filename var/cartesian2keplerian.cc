#include <cmath>
#include <cassert>
#include <cstdio>
#include "ggeodesy/geoconst.hpp"
#include "ggeodesy/units.hpp"

namespace ngpt {
  constexpr double GM = 3.9860044e14;
}

inline void cross_product(const double *a, const double *b, double *c) noexcept {
  c[0] = a[1]*b[2] - a[2]*b[3];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

inline double norm_squared(const double *v) noexcept {
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

inline double dot_product(const double* a, const double *b) noexcept {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
void cartesian2keplerian_(const double *r, const double *rdot) noexcept {
  double hv[3];
  cross_product(r, rdot, hv);
  double h_nsq = norm_squared(hv);
  double hn = std::sqrt(h_nsq);
  double rn = std::sqrt(dot_product(r, r));
  double Omega = std::atan2(hv[0], hv[1]);
  if (Omega<0e0) Omega += ngpt::D2PI;
  assert(Omega>0e0 && Omega<ngpt::D2PI);
  double i = std::acos(hv[2]/hn);
  double rdot2 = dot_product(rdot, rdot);
  rdot2 *= rdot2;
  double Ec = rdot2/2e0 - ngpt::GM / rn;
  double p = h_nsq / ngpt::GM;
  double e = std::sqrt(1e0 + (2e0*h_nsq*Ec)/(ngpt::GM*ngpt::GM));

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
  double h[3];
  cross_product(x, v, h);
  const double h_sq = norm_squared(h);
  const double h_norm = std::sqrt(h_sq);

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
  const double x[] = {5492000.34e0,
                      3984001.40e0,
                         2955.81e0};
  const double v[] = {  -3931.046491e0,
                         5498.676921,
                         3665.980697};

  cartesian2keplerian(x, v);
  cartesian2keplerian_(x, v);
  return 0;
}
