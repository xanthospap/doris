#include "astrodynamics.hpp"
#include <iers2010/matvec.hpp>

void compute_perifocal(double a, double e, double E, double GM, dso::Vector3 &r,
                       dso::Vector3 &v) noexcept {
  const double sE = std::sin(E);
  const double cE = std::cos(E);
  const double me2p1 = std::sqrt((1e0 + e) * (1e0 - e));
  r.x() = a * (cE - e);
  r.y() = a * me2p1 * sE;
  r.z() = 0e0;

  const double vfac = std::sqrt(GM * a) / r.norm();
  v.x() = -vfac * sE;
  v.y() = vfac * me2p1 * cE;
  v.z() = 0e0;
}

int dso::elements2perifocal(const dso::OrbitalElements &ele, double E,
                            double GM, dso::Vector3 &r,
                            dso::Vector3 &v) noexcept {
  const double a = ele.semimajor();
  const double e = ele.eccentricity();

  compute_perifocal(ele.semimajor(), ele.eccentricity(), E, GM, r, v);

  return 0;
}

int dso::elements2perifocal(const dso::OrbitalElements &ele, double GM,
                            dso::Vector3 &r, dso::Vector3 &v) noexcept {
  const double a = ele.semimajor();
  const double e = ele.eccentricity();
  int ok;
  const double E = dso::kepler(e, ele.mean_anomaly(), ok);
  if (ok)
    return ok;
  compute_perifocal(a, e, E, GM, r, v);
  return 0;
}

int dso::perifocal2equatorial(const dso::OrbitalElements &ele,
                              const dso::Vector3 &rp, const dso::Vector3 &vp,
                              dso::Vector3 &re, dso::Vector3 &ve) noexcept {
  const dso::Mat3x3 T =
      perifocal2equatorial_matrix(ele.Omega(), ele.omega(), ele.inclination());
  re = T * rp;
  ve = T * vp;
  return 0;
}

int dso::equatorial2perifocal(const dso::OrbitalElements &ele,
                              const dso::Vector3 &re, const dso::Vector3 &ve,
                              dso::Vector3 &rp, dso::Vector3 &vp) noexcept {
  const dso::Mat3x3 T =
      perifocal2equatorial_matrix(ele.Omega(), ele.omega(), ele.inclination())
          .transpose_inplace();
  rp = T * re;
  vp = T * ve;
  return 0;
}

// transformation is: r_eq = T r_pf
dso::Mat3x3 perifocal2equatorial_matrix(double Omega, double omega, double i) noexcept {
  const double cO = std::cos(Omega);
  const double sO = std::sin(Omega);
  const double co = std::cos(omega);
  const double so = std::sin(omega);
  const double ci = std::cos(i);
  const double si = std::sin(i);


  const double m00 = cO * co - sO * so * ci;
  const double m01 = -cO * so - sO * co * ci;
  const double m02 = sO * si;
  const double m10 = sO * co + cO * so * ci;
  const double m11 = -sO * so + cO * co * ci;
  const double m12 = -cO * si;
  const double m20 = sO * si;
  const double m21 = cO * si;
  const double m22 = ci;
  
  return dso::Mat3x3({m00, m01, m02, m10, m11, m12, m20, m21, m22});
}
