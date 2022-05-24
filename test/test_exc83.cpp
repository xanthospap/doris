#include "astrodynamics.hpp"
#include "datetime/dtcalendar.hpp"
#include "geodesy/car2top.hpp"
#include "geodesy/ellipsoid.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iau.hpp"
#include "iers2010/iersc.hpp"
#include "iers2010/matvec.hpp"
#include <cstdio>

using dso::Vector3;

// Ground station, Bangalore
const Vector3 xsta({1344e3, 6069e3, 1429e3}); // [m]

// zero-epoch for observations
dso::datetime<dso::milliseconds> t0(dso::year(1995), dso::month(3),
                                    dso::day_of_month(30),
                                    dso::milliseconds(0));

// (2-part) state vector, at t0
const Vector3 svr0({-6345e3, -3723e3, -580e3});    // [m]
const Vector3 svv0({2.169e3, -9.266e3, -1.079e3}); // [m/s]

// struct Obs
struct ObsType {
  dso::datetime<dso::milliseconds> t;
  double azim, elev, dist;
};

int main() {

  // artificial observations
  // ObsType obs[6];

  Vector3 r, v;
  // make synthetic/artificial observation set
  for (int i = 0; i < 6; i++) {
    // delta seconds from reference epoch
    const double dt = 1200e0 * (i + 1);

    // epoch ...
    auto t = t0;
    t.add_seconds(dso::milliseconds(static_cast<long>(dt*1e3)));

    // propagate state vector (t0, t)
    // printf("\tinitial    state: %.3f %.3f %.3f %.3f %.3f %.3f\n", svr0.x(),
    // svr0.y(), svr0.z(), svv0.x(), svv0.y(), svv0.z());
    if (dso::propagate_state(iers2010::GMe, svr0, svv0, dt, r, v)) {
      printf("ERROR, failed to propagate state!\n");
      return 1;
    }
    // printf("\tpropagated state: %.3f %.3f %.3f %.3f %.3f %.3f\n", r.x(),
    // r.y(), r.z(), v.x(), v.y(), v.z());

    // site-satellite vector, topocentric (eith limited precision)
    const double gmst = iers2010::sofa::gmst06(t.as_jd(), 0e0, t.as_jd(), 0e0);
    dso::Mat3x3 rz;
    rz.rotz(gmst);
    const Vector3 sat_ec = rz * r;
    Vector3 neu;
    dso::car2top<dso::ellipsoid::grs80>(xsta.x(), xsta.y(), xsta.z(),
                                        sat_ec.x(), sat_ec.y(), sat_ec.z(),
                                        neu.x(), neu.y(), neu.z());
    double s, az, zen;
    dso::top2daz(neu.x(), neu.y(), neu.z(), s, az, zen);

    printf("mjd=%.3f Az=%.3f [deg] Ele=%.3f [deg] R=%.3f [km] GMST=%.6f "
           "[deg] JD=%.5f\n",
           t.as_mjd(), dso::rad2deg(az), dso::rad2deg(dso::DPI - zen), s / 1e3,
           dso::rad2deg(gmst), t.as_jd());
  }

  return 0;
}
