#include "satellite.hpp"
#include <cmath>
#ifdef ASSERT_ERROR
#include <cassert>
#endif
#ifdef DEBUG
#include <cstdio>
#endif

using dso::Vector3;

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
#ifdef DEBUG
  // keplerian[3] = std::sqrt(1e0 + (keplerian[0] / m) * (keplerian[0] / m) *
  //                                   (v * v - (2e0 * m) / r));
  // printf(">> Note difference in magnitudes of e, is: %e\n", keplerian[3]-e);
#endif
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
