#include "astrodynamics.hpp"
#include "iers2010/iersc.hpp"
#include <cmath>
#include <iers2010/matvec.hpp>
#include <limits>
#ifdef ASSERT_ERROR
#include <cassert>
#endif
#ifdef DEBUG
#include <cstdio>
#endif

using dso::Vector3;

int dso::elements2state(const dso::OrbitalElements &ele, double dt, dso::Vector3 &r,
                   dso::Vector3 &v, double GM) noexcept
{
  // mean anomaly M at time t
  const double a = ele.semimajor();
  const double M = ele.mean_anomaly() + std::sqrt(GM/(a*a*a))*dt;

  // solve kepler's equation
  const double e = ele.eccentricity();
  int ok;
  const double E = dso::kepler(e, M, ok);
  // if (ok) return ok;

  // perifocal coordinates (see Montenbruck 2.2.3)
  const double cE = std::cos(E);
  const double sE = std::sin(E);
  const double f = std::sqrt((1e0-e)*(1e0+e)); // aka (1-e^2)^(1/2)
  const double rmag = a*(1e0-e*cE);
  const dso::Vector3 rf({a*(cE-e), a*f*sE, 0e0});
  const double vmag = std::sqrt(GM*a) / rmag;
  const dso::Vector3 vf({-vmag*sE,vmag*f*cE, 0e0});

  // matrix Rz(-Ω)Rx(-i)Rz(-ω)=T , see Gurfil et al, eq. 5.2
  const double sO = std::sin(ele.Omega());
  const double cO = std::cos(ele.Omega());
  const double so = std::sin(ele.omega());
  const double co = std::cos(ele.omega());
  const double si = std::sin(ele.inclination());
  const double ci = std::cos(ele.inclination());
  const dso::Mat3x3 T({cO*co - sO*so*ci, -cO*so - sO*co*ci, sO*si,
                       sO*co + cO*so*ci, -sO*so + cO*co*ci, -cO*si,
                       so*si,             co*si,            ci});
  
  // rotate to get equatorial, aka geocentric/(quasi-)inertial
  r = T * rf;
  v = T * vf;

  return ok;
}

int dso::state2elements(const dso::Vector3 &r, const dso::Vector3 &v,
                   dso::OrbitalElements &ele, double GM) noexcept {
  const Vector3 h = r.cross_product(v);
  const double hmag = h.norm();

  ele.Omega() = std::fmod(std::atan2(h.x(), -h.y()), iers2010::D2PI);
  ele.inclination() = std::atan2(std::sqrt(h.x()*h.x()+h.y()*h.y()), h.y());
  
  // argument of latitude
  const double u = std::atan2(r.z() * hmag, -r.x()*h.y()+r.y()*h.x());
  
  const double rmag = r.norm();
  ele.semimajor() = 1e0 / (2e0/rmag - v.dot_product(v)/GM);
  const double a = ele.semimajor();

  const double eCosE = 1e0 - rmag/a;
  const double eSinE = r.dot_product(v) / std::sqrt(GM*a);
  const double e2 = eCosE*eCosE + eSinE*eSinE;
  
  ele.eccentricity() = std::sqrt(e2);

  // eccentric anomaly
  const double E = std::atan2(eSinE,eCosE);

  ele.mean_anomaly() = std::fmod(E-eSinE, iers2010::D2PI);
  const double true_anomaly = std::atan2(std::sqrt(1e0-e2)*eSinE, eCosE-e2);
  ele.omega() = std::fmod(u-true_anomaly, iers2010::D2PI);

  return 0;
}

int dso::alternatives::state2kepler_montenbruck(const double *state,
                                            double *keplerian,
                                            double m) noexcept {
  // split to position and velocity vectors
  const Vector3 r(Vector3::to_vec3(state));
  const Vector3 v(Vector3::to_vec3(state + 3));
  const double rmag = r.norm();

  const Vector3 h = r.cross_product(v);
  const double hmag = h.norm();
  keplerian[0] = hmag;
  const Vector3 W{{h.x()/hmag, h.y()/hmag, h.z()/hmag}};

  keplerian[1] = std::atan2(std::sqrt(W.x()*W.x()+W.y()*W.y()), W.z());
  keplerian[2] = std::atan2(W.x(), -W.y());

  const double p = h.norm_squared() / m;
  const double a = 1e0 / ( (2e0/rmag) - (v.norm_squared() / m) );
  const double n = std::sqrt(m / (a*a*a));
  const double e = std::sqrt(1e0 - (p/a));
  keplerian[3] = e;

  const double E = atan2(r.dot_product(v)/a*a/n, 1e0-rmag/a);
  const double sE = std::sin(E);
  // const double M = E - e * sE;

  const double u = std::atan2(r.z(), -r.x()*W.y() + r.y()*W.x());
  const double t = std::atan2(std::sqrt(1e0-e*e) * sE, std::cos(E)-e);
  keplerian[4] = u - t;
  keplerian[5] = t;

  return 0;
}

int dso::alternatives::state2kepler_vallado(const double *state,
                                            double *keplerian,
                                            double m) noexcept {
  // split to position and velocity vectors
  const Vector3 r(Vector3::to_vec3(state));
  const Vector3 v(Vector3::to_vec3(state + 3));
  const double rnorm = r.norm();

  const Vector3 h = r.cross_product(v);
  const double hnorm = h.norm();
  keplerian[0] = hnorm;

  const Vector3 n{{-h.y(), h.x(), 0e0}}; // N = k_unit x h
  const double nnorm = n.norm();

  Vector3 e = r * (v.norm_squared() - m / rnorm) - v * (r.dot_product(v));
  e /= m;
  double enorm = e.norm();
  keplerian[3] = enorm;

//  const double ksi = v.norm_squared() / 2e0 - m / rnorm;
//#ifdef BRANCHLESS
//  double alpha = -m / (1e0 * ksi) +
//          (enorm == 1e0) * std::numeric_limits<double>::infinity(); // [return]
//  double rho = (enorm != 1e0) * (alpha * (1e0 - enorm * enorm)) +
//        (enorm == 1e0) * (h.norm_squared() / m);
//#else
//  double alpha =
//      (enorm == 1e0) ? std::numeric_limits<double>::infinity() : -m / (1e0 * ksi);
//  double rho = (enorm == 1e0) ? (h.norm_squared() / m) : alpha * (1e0 - enorm * enorm);
//#endif

  keplerian[1] = std::acos(h.z() / hnorm);

  double Omega = std::acos(n.x() / nnorm);
  if (n.y() < 0e0)
    Omega = iers2010::D2PI - Omega;
  keplerian[2] = Omega;

  double omega = std::acos(n.dot_product(e) / (nnorm * enorm));
  if (e.z() < 0e0)
    omega = iers2010::D2PI - omega;
  keplerian[4] = omega;

  double ni = std::acos(e.dot_product(r) / (enorm * rnorm));
  if (r.dot_product(v) < 0e0)
    ni = iers2010::D2PI - ni;
  keplerian[5] = ni;

  return 0;
}

int dso::state2kepler(const double *state, double *keplerian,
                      double m) noexcept {
  // split to position and velocity vectors
  Vector3 pos(Vector3::to_vec3(state));
  Vector3 vel(Vector3::to_vec3(state + 3));

  // calculate distance and speed
  const double r = pos.norm();
  const double v = vel.norm();

  // calculate radial velocity, v_r = (r * v) / |r|
  // Note that if vr > 0, the spacecraft is flying away from perigee. If
  // vr < 0, it is flying toward perigee.
  const double vr = pos.dot_product(vel) / r;

  // Calculate the specific angular momentum, h = r x v
  const Vector3 hv = pos.cross_product(vel);
  // and its magnitude (h); This is the first orbital element.
  keplerian[0] = hv.norm();

  // Calculate the inclination (i); this is the second orbital element.
  // Recall that i must lie between 0° and 180°, which is precisely the range
  // (principal values) of the acos function
  // (https://en.cppreference.com/w/cpp/numeric/math/acos). Hence, there is no
  // quadrant ambiguity to contend with here. If 90° < i  180°, the angular
  // momentum h points in a southerly direction. In that case, the orbit is
  // retrograde, which means that the motion of the satellite around the earth
  // is opposite to earth’s rotation.
  keplerian[1] = std::acos(hv.z() / keplerian[0]);

  // calculate N = k x h ; this vector defines the node line.
  const Vector3 Nv{{-hv.y(), hv.x(), 0e0}};
  const double N = Nv.norm();

  // Calculate the right ascension of the ascending node (Omega).
  // This is the third orbital element. If NX>0, then Ω lies in either the
  // first or fourth quadrant. If NX < 0, then Ω lies in either the second or
  // third quadrant. To place Ω in the proper quadrant, observe that the
  // ascending node lies on the positive side of the vertical XZ plane
  // (0 <  Ω < 180°) if NY > 0. On the other hand, the ascending node lies on
  // the negative side of the XZ plane (180° < Ω < 360°) if NY < 0. Therefore,
  // NY > 0 implies that 0  Ω < 180°, whereas NY < 0 implies that
  // 180° < Ω < 360°
  const double Nf = std::acos(Nv.x() / N);
#ifdef BRANCHLESS
  const double retN[] = {Nf, iers2010::D2PI - Nf};
  keplerian[2] = retN[Nv.y() < 0e0];
#else
  keplerian[2] = (Nv.y() < 0e0) ? (iers2010::D2PI - Nf) : Nf;
#endif

  // Calculate the eccentricity vector e
  const Vector3 ev = (pos * (v * v - (m / r)) - vel * (r * vr)) / m;
  // or ... just compute the magnitude ...., the forth orbital element
  // keplerian[3] = std::sqrt(1e0 + (h/m)*(h/m)*(v*v - (2e0*m) / r));
  const double e = ev.norm();
  keplerian[3] = e;

  // Calculate the argument of perigee.
  // This is the fifth orbital element. If Ne > 0, then ω lies in either the
  // first or fourth quadrant. If Ne < 0, then ω lies in either the second or
  // third quadrant. To place ω in the proper quadrant, observe that perigee
  // lies above the equatorial plane (0°  ω < 180°) if e points up (in the
  // positive Z direction) and that perigee lies below the plane
  // (180°  ω < 360°) if e points down. Therefore, eZ >=  0 implies that
  // 0° < ω < 180°, whereas eZ < 0 implies that 180° < ω < 360°.
  const double Of = std::acos(Nv.dot_product(ev) / N / e);
#ifdef BRANCHLESS
  const double retO[] = {Of, iers2010::D2PI - Of};
  keplerian[4] = retO[ev.z() < 0e0];
#else
  keplerian[4] = (ev.z() < 0e0) ? (iers2010::D2PI - Of) : Of;
#endif

  // Calculate true anomaly.
  // This is the sixth and final orbital element. If e*r > 0, then θ lies in
  // the first or fourth quadrant. If e*r < 0, then θ lies in the second or
  // third quadrant. To place θ in the proper quadrant, note that if the
  // satellite is flying away from perigee (r*v >= 0), then 0 <  θ < 180°,
  // whereas if the satellite is flying toward perigee (r*v < 0), then
  // 180° < θ < 360°.
  const double Tf = std::acos(ev.dot_product(pos) / e / r);
#ifdef BRANCHLESS
  const double retT[] = {Tf, iers2010::D2PI - Tf};
  keplerian[5] = retT[vr < 0e0];
#else
  keplerian[5] = (vr < 0e0) ? (iers2010::D2PI - Tf) : Tf;
#endif

  return 0;
}
