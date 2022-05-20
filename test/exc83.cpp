#include <cstdio>
#include "datetime/dtcalendar.hpp"
#include "iers2010/iersc.hpp"
#include "iers2010/matvec.hpp"
#include "astrodynamics.hpp"

using dso::Vector3;

// Ground station, Bangalore
const Vector3 xsta({1344e3,6069e3,1429e3}); // [m]

// kinda zero-epoch for observations
dso::datetime<dso::milliseconds> t0(ngpt::year(1995), ngpt::month(3),
                                    ngpt::day_of_month(30),
                                    dso::milliseconds(0));

// (2-part) state vector, at t0
const Vector3 svr0({-6345e3, -3723e3, -580e3}); // [m]
const Vector3 svv0({2.169e3, -9.266e3, -1.079e3}); // [m/s]

// struct Obs
struct ObsType {
  dso::datetime<dso::milliseconds> t;
  double azim, elev, dist;
};

int main() {

  // artificial observations
  ObsType obs[6];

  Vector3 r,v;
  // make synthetic/artificial observation set
  for (int i=0; i<6; i++) {
    // epoch ...
    auto t = t0;
    t.add_seconds(dso::milliseconds(1'200'000L*(i+1)));

    // delta seconds from reference epoch
    const double dt = 1200e0 * (i+1);

    // propagate state vector (t0, t)
    int error = dso::propagate_state(iers2010::GMe, svr0, svv0, dt, r, v);
    
  }
}
