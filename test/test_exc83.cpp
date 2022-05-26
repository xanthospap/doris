#include "astrodynamics.hpp"
#include "datetime/datetime_write.hpp"
#include "datetime/dtcalendar.hpp"
#include "geodesy/car2top.hpp"
#include "geodesy/ellipsoid.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iau.hpp"
#include "iers2010/iersc.hpp"
#include "iers2010/matvec.hpp"
#include "filters/ekf.hpp"
#include <cstdio>
#include <datetime/dtfund.hpp>

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

  // a buffer to write datetime info
  char buf[64];

  // artificial observations
  ObsType obs[6];

  Vector3 r, v;
  printf("Date       A [deg] Ele [deg] Range [km]\n");
  // make synthetic/artificial observation set
  for (int i = 0; i < 6; i++) {
    // delta seconds from reference epoch
    const double dt = 1200e0 * (i + 1);

    // epoch ...
    auto t = t0;
    t.add_seconds(dso::milliseconds(static_cast<long>(dt * 1e3)));

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

    printf("%s %.3f %.3f %.3f\n", dso::strftime_ymd_hmfs(t, buf),
           dso::rad2deg(az), dso::rad2deg(dso::DPI/2e0 - zen), s / 1e3);

    obs[i] = ObsType({t, az, dso::DPI/2e0 - zen, s});
  }

  // Orbit determination; initial time t0
  constexpr const int N = 6;
  Eigen::Matrix<double, N, 1> Y;
  Y(0) = svr0.x() + 10e3; Y(1) = svr0.y() - 5e3; Y(2) = svr0.z()+1e3;
  Y(3) = svv0.x() - 1e0; Y(4) = svr0.y() + 3e0; Y(5) = svr0.z() - 5e-1;
  
  Eigen::Matrix<double, N, N> P = Eigen::Matrix<double, N, N>::Zeros();
  P(0,0) = P(1,1) = P(2,2) = 1e8;
  P(3,3) = P(4,4) = P(5,5) = 1e2;

  dso::ExtendedKalmanFilter<6,dso::milliseconds> Filter(t0, Y, P);

  // measurements ...
  for (int i=0; i<6; i++) {

    // previous step
    auto t_prev = Filter.t;
    auto state_prev = Filter.x;

    // propagation to measurement epoch
    auto t = Obs[i].t;
    Vector3 yr({state_prev(0), state_prev(1), state_prev(2)});
    Vector3 yv({state_prev(3), state_prev(4), state_prev(5)});
    double dt = (t.delta_sec(t0)).to_fractional_seconds();
    dso::propagate_state(iers2010::GMe, yr, yv, dt, r, v);
    auto Phi = dso::state_partials(iers2010::GMe,

    Y(0) = r.x(); Y(1) = r.y(); Y(2) = r.z();
    Y(3) = v.x(); Y(4) = v.y(); Y(5) = v.z();
  }

  return 0;
}
