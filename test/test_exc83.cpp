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

// struct Obs
struct ObsType {
  dso::datetime<dso::milliseconds> t;
  double azim, elev, dist;
};

constexpr const double sigma_range = 10e0;               // [m]
constexpr const double sigma_angle = dso::deg2rad(1e-2); // [rad]
constexpr const double GM = /*iers2010::GMe;*/ 398600.4415e+9;

int main() {
  // zero-epoch for observations
  dso::datetime<dso::milliseconds> t0(dso::year(1995), dso::month(3),
                                      dso::day_of_month(30),
                                      dso::milliseconds(0));

  // state vector, at t0 [m] and [m/s]
  Eigen::Matrix<double, 6, 1> Y0, Y;
  Y0 << -6345e3, -3723e3, -580e3, 2.169e3, -9.266e3, -1.079e3;
  // Ground station, Bangalore
  Eigen::Matrix<double, 3, 1> xsta;
  xsta << 1344e3, 6069e3, 1429e3; //[m]

  // a buffer to write datetime info
  char buf[64];

  // artificial observations
  ObsType obs[6];

  printf("Date       A [deg] Ele [deg] Range [km]\n");
  // make synthetic/artificial observation set
  for (int i = 0; i < 6; i++) {
    // delta seconds from reference epoch
    const double dt = 1200e0 * (i + 1);

    // epoch ...
    auto t = t0;
    t.add_seconds(dso::milliseconds(static_cast<long>(dt * 1e3)));

    Y = dso::propagate_state(GM, Y0, dt);
    // printf("\tpropagated state to [%.3f %.3f %.3f %.5f %.5f
    // %.5f]\n",Y(0),Y(1),Y(2),Y(3),Y(4),Y(5));

    // site-satellite vector, topocentric (eith limited precision)
    const double gmst = iers2010::sofa::gmst06(t.as_jd(), 0e0, t.as_jd(), 0e0);

    // Eigen::Matrix<double, 3, 3> Rz(Eigen::Matrix<double, 3, 3>::Identity());
    Eigen::Matrix<double, 3, 3> Rz(
        Eigen::AngleAxisd(-gmst, Eigen::Vector3d::UnitZ()));
    Eigen::Matrix<double, 3, 1> enu =
        dso::car2top<dso::ellipsoid::grs80>(xsta, Rz * Y.head<3>());
    // printf("\tMatrix Rz= | %.3f %.3f %.3f |\n", Rz(0,0),Rz(0,1),Rz(0,2));
    // printf("\t           | %.3f %.3f %.3f |\n", Rz(1,0),Rz(1,1),Rz(1,2));
    // printf("\t           | %.3f %.3f %.3f |\n", Rz(2,0),Rz(2,1),Rz(2,2));
    // auto RzY = Rz * Y.head<3>();
    // printf("\tXYZ = | %.3f %.3f %.3f | (sta)\n", xsta(0),xsta(1),xsta(2));
    // printf("\tXYZ = | %.3f %.3f %.3f | (sat)\n", RzY(0),RzY(1),RzY(2));
    // printf("\tENU = | %.3f %.3f %.3f |\n", enu(0),enu(1),enu(2));

    double s, az, el;
    s = dso::top2dae(enu, az, el);

    printf("%s %.3f %.3f %.3f\n", dso::strftime_ymd_hmfs(t, buf),
           dso::rad2deg(az), dso::rad2deg(el), s / 1e3);

    obs[i] = ObsType({t, az, el, s});
  }

  // Orbit determination; initial time t0
  constexpr const int N = 6;
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
  //printf("LTC Matrix:\n");
  //for (int i = 0; i < 3; i++) {
  //  printf("  |");
  //  for (int k = 0; k < 3; k++)
  //    printf("%+.3e ", lct(i, k));
  //  printf("|\n");
  //}

  // measurements ...
  double s, az, el;
  Eigen::Matrix<double, 3, 1> dAdenu, dEdenu;
  for (int i = 0; i < 6; i++) {

    // printf("Orbit determination, iteration %d\n", i);

    // previous step
    auto Yprev = Filter.state();
    auto tprev = Filter.time();

    // propagation to measurement epoch
    auto t = obs[i].t;
    double dt = (t.delta_sec(tprev)).to_fractional_seconds();
    //printf(
    //    "\tPropagating state = [%.1f %.1f %.1f %.3f %.3f %.3f] with dt=%.1f\n",
    //    Yprev(0), Yprev(1), Yprev(2), Yprev(3), Yprev(4), Yprev(5), dt);
    Y = dso::propagate_state(GM, Yprev, dt, Phi);
    //printf("\tPropagated  state = [%.1f %.1f %.1f %.3f %.3f %.3f]\n", Y(0),
    //       Y(1), Y(2), Y(3), Y(4), Y(5));
    // printf("\tpropagation to measurement epoch:");
    // printf("\tdt=%.5f Y=[%.3f %.3f %.3f %.5f %.5f %.5f]\n",
    // dt,Y(0),Y(1),Y(2),Y(3),Y(4),Y(5));

    // Earth rotation angle to be used for local tangent coords
    const double gmst = iers2010::sofa::gmst06(t.as_jd(), 0e0, t.as_jd(), 0e0);
    Eigen::Matrix<double, 3, 3> Rz(Eigen::Matrix<double, 3, 3>::Identity());
    Rz = Eigen::AngleAxisd(-gmst, Eigen::Vector3d::UnitZ());

    // time update
    Filter.time_update(t, Y, Phi);

    // truth orbit
    Ytrue = dso::propagate_state(GM, Y0, dt, Phi_true);

    // Handle azimouth measurement
    vec3 enu = dso::car2top<dso::ellipsoid::grs80>(
        xsta, Rz * Filter.state().head<3>());
    s = dso::top2dae(enu, az, el, dAdenu, dEdenu);
    // printf("\tdAdenu=[%.3e %.3e %.3e]\n", dAdenu(0), dAdenu(1), dAdenu(2));
    Eigen::Matrix<double, 6, 1> dAdY;
    dAdY << (dAdenu.transpose() * lct * Rz).transpose(),
        Eigen::Matrix<double, 3, 1>::Zero();
    // printf("\tdAdY  =[%.3e %.3e %.3e]\n", dAdY(0), dAdY(1), dAdY(2));

    // Measurement update (azimouth)
    Filter.observation_update(obs[i].azim, az, sigma_angle / std::cos(el),
                              dAdY);
    //printf("\tUpdated state     = [%.1f %.1f %.1f %.3f %.3f %.3f] (azimouth)\n",
    //       Filter.state()(0), Filter.state()(1), Filter.state()(2),
    //       Filter.state()(3), Filter.state()(4), Filter.state()(5));

    // Handle elevation measurement
    enu = dso::car2top<dso::ellipsoid::grs80>(xsta,
                                              Rz * Filter.state().head<3>());
    s = dso::top2dae(enu, az, el, dAdenu, dEdenu);
    // printf("\tdEdenu=[%.3e %.3e %.3e]\n", dEdenu(0), dEdenu(1), dEdenu(2));
    Eigen::Matrix<double, 6, 1> dEdY;
    dEdY << (dEdenu.transpose() * lct * Rz).transpose(),
        Eigen::Matrix<double, 3, 1>::Zero();

    // Measurement update (elevation)
    Filter.observation_update(obs[i].elev, el, sigma_angle, dEdY);
    //printf(
    //    "\tUpdated state     = [%.1f %.1f %.1f %.3f %.3f %.3f] (elevation)\n",
    //    Filter.state()(0), Filter.state()(1), Filter.state()(2),
    //    Filter.state()(3), Filter.state()(4), Filter.state()(5));

    // Handle range measurement
    enu = dso::car2top<dso::ellipsoid::grs80>(xsta,
                                              Rz * Filter.state().head<3>());
    s = dso::top2dae(enu, az, el);
    Eigen::Matrix<double, 3, 1> dDdenu = enu / enu.norm();
    //printf("\tdDdenu=[%.3e %.3e %.3e]\n", dDdenu(0), dDdenu(1), dDdenu(2));
    Eigen::Matrix<double, 6, 1> dDdY;
    dDdY << (dDdenu.transpose() * lct * Rz).transpose(),
        Eigen::Matrix<double, 3, 1>::Zero();
    //printf("\tdDdY  =[%.3e %.3e %.3e]\n", dDdY(0), dDdY(1), dDdY(2));

    // Measurement update (elevation)
    Filter.observation_update(obs[i].dist, s, sigma_range, dDdY);
    //printf("\tUpdated state     = [%.1f %.1f %.1f %.3f %.3f %.3f] (range)\n",
    //       Filter.state()(0), Filter.state()(1), Filter.state()(2),
    //       Filter.state()(3), Filter.state()(4), Filter.state()(5));

    // auto x = Filter.state();
    // printf("State: [%.3f %.3f %.3f] m.\n", x(0), x(1), x(2));
    // printf("       [%.5f %.5f %.5f] m/s.\n", x(3), x(4), x(5));
  }
  printf("Final Result:\n");
  for (int i = 0; i < 6; i++)
    printf("\t%+20.6f +- %15.6f\n", Filter.state()(i),
           std::sqrt(Filter.P(i, i)));

  return 0;
}
