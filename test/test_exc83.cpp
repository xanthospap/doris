#include "astrodynamics.hpp"
#include "datetime/datetime_write.hpp"
#include "datetime/dtcalendar.hpp"
#include "filters/ekf.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iau.hpp"
#include "iers2010/iersc.hpp"
#include "matvec/matvec.hpp"
#include <cstdio>

// using dso::Vector3;
using vec3 = Eigen::Matrix<double, 3, 1>;
using vec6 = Eigen::Matrix<double, 6, 1>;

// Eigen::Matrix<double,6,1> rv2state(const Vector3 &r, const Vector3 &v)
// noexcept {
//   Eigen::Matrix<double,6,1> Y;
//   Y(0,0) = r.x(); Y(0,1) = r.y(); Y(0,2) = r.z();
//   Y(0,3) = v.x(); Y(0,4) = v.y(); Y(0,5) = v.z();
//   return Y;
// }

// Ground station, Bangalore
// const Vector3 xsta({1344e3, 6069e3, 1429e3}); // [m]

// state vector, at t0 [m] and [m/s]
// vec6 Y0; Y0 << -6345e3, -3723e3, -580e3, 2.169e3, -9.266e3, -1.079e3;
// const Vector3 svr0({-6345e3, -3723e3, -580e3});    // [m]
// const Vector3 svv0({2.169e3, -9.266e3, -1.079e3}); // [m/s]

// struct Obs
struct ObsType {
  dso::datetime<dso::milliseconds> t;
  double azim, elev, dist;
};

constexpr const double sigma_range = 10e0; // [m]
constexpr const double sigma_angle = dso::rad2deg(1e-2); // [rad]
constexpr const double GM = iers2010::GMe;

int main() {
  // zero-epoch for observations
  dso::datetime<dso::milliseconds> t0(dso::year(1995), dso::month(3),
                                      dso::day_of_month(30),
                                      dso::milliseconds(0));

  // Ground station, Bangalore
  Eigen::Matrix<double, 3, 1> xsta;
  // state vector, at t0 [m] and [m/s]
  vec6 Y0;
  Y0 << -6345e3, -3723e3, -580e3, 2.169e3, -9.266e3, -1.079e3;
  xsta << 1344e3, 6069e3, 1429e3; //[m]

  // a buffer to write datetime info
  char buf[64];

  // artificial observations
  ObsType obs[6];

  vec3 r, v;
  vec6 Y;
  printf("Date       A [deg] Ele [deg] Range [km]\n");
  // make synthetic/artificial observation set
  for (int i = 0; i < 6; i++) {
    // delta seconds from reference epoch
    const double dt = 1200e0 * (i + 1);

    // epoch ...
    auto t = t0;
    t.add_seconds(dso::milliseconds(static_cast<long>(dt * 1e3)));

    Y = dso::propagate_state(iers2010::GMe, Y0, dt);

    // site-satellite vector, topocentric (eith limited precision)
    const double gmst = iers2010::sofa::gmst06(t.as_jd(), 0e0, t.as_jd(), 0e0);

    // dso::Mat3x3 rz;
    // rz.rotz(gmst);
    // const Vector3 sat_ec = rz * r;
    // printf("Note enu=[%.3f,%.3f,%.3f]\n", enu(0),enu(1),enu(2));

    Eigen::Matrix<double, 3, 3> Rz(Eigen::Matrix<double, 3, 3>::Identity());
    Rz = Eigen::AngleAxisd(gmst, Eigen::Vector3d::UnitZ());
    vec3 enu = dso::car2top<dso::ellipsoid::grs80>(
        xsta, /*sat_ec*/ Rz.transpose() * Y.head<3>());
    // auto foo = Rz.transpose() * Y.head<3>();
    // printf("\tcar2top for (%.3f %.3f %.3f) - (%.3f %.3f %.3f)\n", foo(0),
    // foo(1), foo(2), xsta(0), xsta(1), xsta(2)); printf("\tENU=(%.3f %.3f
    // %.3f)\n", enu(0),enu(1),enu(2));

    double s, az, el;
    dso::top2dae(enu, s, az, el);

    // printf("\tState=(%.3f %.3f %.3f %.3f %.3f %.3f) [m,m/s]\n\tgmst=%.5f
    // deg.\n", Y(0),Y(1),Y(2),Y(3),Y(4),Y(5),dso::rad2deg(gmst));
    printf("%s %.3f %.3f %.3f\n", dso::strftime_ymd_hmfs(t, buf),
           dso::rad2deg(az), dso::rad2deg(el), s / 1e3);

    obs[i] = ObsType({t, az, el, s});
  }

  // Orbit determination; initial time t0
  constexpr const int N = 6;
  // Eigen::Matrix<double, N, 1> Y;
  Y = Y0;
  Y(0) += 10e3;
  Y(1) += -5e3;
  Y(2) += 1e3;
  Y(3) += -1e0;
  Y(4) += 3e0;
  Y(5) += -5e-1;

  Eigen::Matrix<double, N, N> P = Eigen::Matrix<double, N, N>::Zero();
  P(0, 0) = P(1, 1) = P(2, 2) = 1e8;
  P(3, 3) = P(4, 4) = P(5, 5) = 1e2;

  dso::ExtendedKalmanFilter<6, dso::milliseconds> Filter(t0, Y, P);
  Eigen::Matrix<double, 6, 6> Phi;

  Eigen::Matrix<double, 6, 6> Phi_true;
  Eigen::Matrix<double, N, 1> Ytrue;

  Eigen::Matrix<double, 3, 3> lct(
      dso::topocentric_matrix(dso::car2ell<dso::ellipsoid::grs80>(xsta)));

  // measurements ...
  double s, az, el;
  Eigen::Matrix<double, 3, 1> dAdenu, dEdenu;
  for (int i = 0; i < 6; i++) {

    // previous step
    // auto t_prev = Filter.t;
    auto Yprev = Filter.x;

    // propagation to measurement epoch
    auto t = obs[i].t;
    // Vector3 yr({state_prev(0), state_prev(1), state_prev(2)});
    // Vector3 yv({state_prev(3), state_prev(4), state_prev(5)});
    double dt = (t.delta_sec(t0)).to_fractional_seconds();
    // dso::propagate_state(GM, yr, yv, dt, r, v, Phi);
    Y = dso::propagate_state(GM, Yprev, dt, Phi);

    // time update
    Filter.time_update(t, Y, Phi);

    // truth orbit
    Ytrue = dso::propagate_state(GM, Y0, dt, Phi_true);

    // Handle azimouth measurement
    const double gmst = iers2010::sofa::gmst06(t.as_jd(), 0e0, t.as_jd(), 0e0);
    Eigen::Matrix<double, 3, 3> Rz(Eigen::Matrix<double, 3, 3>::Identity());
    Rz = Eigen::AngleAxisd(gmst, Eigen::Vector3d::UnitZ());
    vec3 enu =
        dso::car2top<dso::ellipsoid::grs80>(xsta, Rz.transpose() * Y.head<3>());
    dso::top2dae(enu, s, az, el, dAdenu, dEdenu);
    Eigen::Matrix<double, 6, 1> dAdY;
    dAdY << (dAdenu.transpose() * lct * Rz).transpose(),
        Eigen::Matrix<double, 3, 1>::Zero();

    // Measurement update (azimouth)
    Filter.observation_update(obs[i].azim,az,sigma_angle/std::cos(el),dAdY);
  }

  return 0;
}
