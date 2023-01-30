#include "astrodynamics.hpp"

// transformation is: r_eq = T r_pf
Eigen::Matrix<double,3,3> dso::perifocal2equatorial_matrix(double Omega, double omega,
                                             double i) noexcept {
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

  Eigen::Matrix<double,3,3> R;
  R << m00, m01, m02, m10, m11, m12, m20, m21, m22;
  return R;
}

/*
void compute_perifocal(double a, double e, double E, double GM, dso::Vector3 &r,
                       dso::Vector3 &v) noexcept {
  const double me2p1 = std::sqrt((1e0 + e) * (1e0 - e));
  r.x() = a * (cE - e);
  r.y() = a * me2p1 * sE;
  r.z() = 0e0;

  const double vfac = std::sqrt(GM * a) / r.norm();
  v.x() = -vfac * sE;
  v.y() = vfac * me2p1 * cE;
  v.z() = 0e0;
}*/

/*
int dso::elements2perifocal(const dso::OrbitalElements &ele, double E,
                            double GM, dso::Vector3 &r,
                            dso::Vector3 &v) noexcept {
  compute_perifocal(ele.semimajor(), ele.eccentricity(), E, GM, r, v);
  return 0;
}*/
/*
Eigen::Matrix<double, 6, 1>
dso::elements2perifocal(const dso::OrbitalElements &ele, double E,
                        double GM) noexcept {
  return compute_perifocal(ele.semimajor(), ele.eccentricity(), E, GM);
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
*/
Eigen::Matrix<double, 6, 1> dso::core::elements2perifocal(double GM, double E,
                                                          double e,
                                                          double a) noexcept {
  const double sE = std::sin(E);
  const double cE = std::cos(E);
  const double me2p1 = std::sqrt((1e0 + e) * (1e0 - e));
  const double r = a * (1e0 - e * cE); // Distance
  const double v = std::sqrt(GM * a) / r;
  Eigen::Matrix<double, 6, 1> P;
  P << a * (cE - e), a * me2p1 * sE, 0e0, -v * sE, v * me2p1 * cE, 0e0;
  return P;
}
Eigen::Matrix<double, 6, 1>
dso::elements2perifocal(double GM, const dso::OrbitalElements &ele,
                        double dt) noexcept {
  const double a = ele.semimajor();
  const double e = ele.eccentricity();
  const double n = std::sqrt(GM / (a * a * a));
  const double M = ele.mean_anomaly() + n * dt;
  int error;
  const double E = dso::kepler(e, M, error);
  return core::elements2perifocal(GM, E, e, a);
}
Eigen::Matrix<double, 6, 1>
dso::elements2perifocal(double GM, const dso::OrbitalElements &ele) noexcept {
  const double a = ele.semimajor();
  const double e = ele.eccentricity();
  int error;
  const double E = dso::kepler(e, ele.mean_anomaly(), error);
  return core::elements2perifocal(GM, E, e, a);
}

int dso::perifocal2equatorial(const dso::OrbitalElements &ele,
                              const Eigen::Matrix<double,3,1> &rp, const Eigen::Matrix<double,3,1> &vp,
                              Eigen::Matrix<double,3,1> &re, Eigen::Matrix<double,3,1> &ve) noexcept {
  const Eigen::Matrix<double,3,3> T =
      perifocal2equatorial_matrix(ele.Omega(), ele.omega(), ele.inclination());
  re = T * rp;
  ve = T * vp;
  return 0;
}
Eigen::Matrix<double, 6, 1>
dso::perifocal2equatorial(const dso::OrbitalElements &ele,
                          const Eigen::Matrix<double, 6, 1> &Yperif) noexcept {
  Eigen::Matrix<double,3,1> rp;rp<<Yperif(0), Yperif(1), Yperif(2);
  Eigen::Matrix<double,3,1> vp;vp<<Yperif(3), Yperif(4), Yperif(5);
  Eigen::Matrix<double,3,1> r, v;
  dso::perifocal2equatorial(ele, rp, vp, r, v);
  Eigen::Matrix<double, 6, 1> Y;
  Y << r(0), r(1), r(2), v(0), v(1), v(2);
  return Y;
}

int dso::equatorial2perifocal(const dso::OrbitalElements &ele,
                              const Eigen::Matrix<double,3,1> &re, const Eigen::Matrix<double,3,1> &ve,
                              Eigen::Matrix<double,3,1> &rp, Eigen::Matrix<double,3,1> &vp) noexcept {
  const Eigen::Matrix<double,3,3> T =
      perifocal2equatorial_matrix(ele.Omega(), ele.omega(), ele.inclination())
          .transpose();
  rp = T * re;
  vp = T * ve;
  return 0;
}
Eigen::Matrix<double, 6, 1>
dso::equatorial2perifocal(const dso::OrbitalElements &ele,
                          const Eigen::Matrix<double, 6, 1> &Yeq) noexcept {
  Eigen::Matrix<double,3,1> req;req<<Yeq(0), Yeq(1), Yeq(2);
  Eigen::Matrix<double,3,1> veq;veq<<Yeq(3), Yeq(4), Yeq(5);
  Eigen::Matrix<double,3,1> r, v;
  dso::equatorial2perifocal(ele, req, veq, r, v);
  Eigen::Matrix<double, 6, 1> Y;
  Y << r(0), r(1), r(2), v(0), v(1), v(2);
  return Y;
}
