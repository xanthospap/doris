#include "astrodynamics.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "orbit_integration.hpp"
#include "iers2010/cel2ter.hpp"

constexpr const int m = 0;
constexpr const int n = 6 + m;
/*
 * tsec : seconds from reference time
 *  aka current time t = t0 + tsec
 * 
 * yP0 : Initial condition for state and state transition matrix, in column-
 *  wise oder, i.e.
 *  yP0 = [ y0, Phi(t0,t0) ] 
 *      = [ r0_(3x1), v0_(3x1), 
 *          dr/dr0_(3x3), dr/dv0_(3x3), 
 *          dv/dr0_(3x3), dv/dv0_(3x3) 
 *        ]
 *      = [ x0, y0, z0, v_x0, v_y0, v_z0,             : 6 (elements)
 *          { dx/dx0, dx/dy0, dx/dz0,
 *            dy/dx0, dy/dy0, dy/dz0,
 *            dz/dx0, dz/dy0, dz/dz0 },               : +9 [dr / dr0]
 *          { dx/dv_x0, dx/dv_y0, dx/dv_z0,
 *            dy/dv_x0, dy/dv_y0, dy/dv_z0,
 *            dz/dv_x0, dz/dv_y0, dz/dv_z0 },         : +9 [dr / dv0]
 *          { dv_x/dx0, dv_x/dy0, dv_x/dz0,
 *            dv_y/dx0, dv_y/dy0, dv_y/dz0,
 *            dv_z/dx0, dv_z/dy0, dv_z/dz0 },         : +9 [dv / dr0]
 *          { dv_x/dv_x0, dv_x/dv_y0, dv_x/dv_z0,
 *            dv_y/dv_x0, dv_y/dv_y0, dv_y/dv_z0,
 *            dv_z/dv_x0, dv_z/dv_y0, dv_z/dv_z0 }    : +9 [dv / dv0]
 *        ] aka a vector of 42 elements/rows
 *  Initially, the state transition matrix should be set to the identity
 *  matrix, i.e. Phi(t0,t0) = I_(6x6)
 */
void dso::VariationalEquations(
    // seconds from reference epoch (TAI)
    double tsec,
    // state and state transition matrix (inertial RF)
    const Eigen::VectorXd &yP0,
    // state derivative and state transition matrix derivative (inertial RF)
    Eigen::Ref<Eigen::VectorXd> yPt,
    // auxiliary parametrs
    dso::IntegrationParameters &params) noexcept {

  // current mjd, TAI
  const dso::TwoPartDate _cmjd(params.mjd_tai +
                              dso::TwoPartDate(0e0, tsec / 86400e0));
  const dso::TwoPartDate cmjd = _cmjd.normalized();

  // terrestrial to celestial for epoch
  dso::Itrs2Gcrs Rot(cmjd.tai2tt(), &params.eopLUT);

  // split position and velocity vectors (inertial)
  Eigen::Matrix<double, 3, 1> r = yP0.block<3, 1>(0, 0);
  Eigen::Matrix<double, 3, 1> v = yP0.block<3, 1>(3, 0);

  // f = y''(t_0,y_0) = a(t_0, y_0) [=(a_x, a_y, a_z)]
  Eigen::Matrix<double, 3, 1> f  = Eigen::Matrix<double, 3, 1>::Zero();
  // df = grav(f(t_0,y_0)) = grad(a(t_0, y_0))  
  //      | da_x/dx0 da_y/dy0 da_z/dz0 |
  //    = | da_x/dx0 da_y/dy0 da_z/dz0 |
  //      | da_x/dx0 da_y/dy0 da_z/dz0 |
  Eigen::Matrix<double, 3, 3> df = Eigen::Matrix<double, 3, 3>::Zero();

  // satellite position at t, in ECEF
  Eigen::Matrix<double, 3, 1> r_geo = Rot.gcrf2itrf(r);
  { // compute gravity-induced acceleration (we need the position vector in
    // ITRF)
    Eigen::Matrix<double, 3, 3> gpartials;
    Eigen::Matrix<double, 3, 1> gacc;
    dso::gravity_acceleration(params.harmonics, r_geo, params.degree,
                   params.harmonics.Re(), params.harmonics.GM(), gacc,
                   gpartials, params.V, params.W);

    // gravity acceleration in earth-fixed frame; need to have inertial 
    // acceleration!
    gacc = Rot.itrf2gcrf(gacc);
    const auto T = Rot.itrf2gcrf();
    gpartials = T * gpartials * T.transpose();

    f += gacc;
    df += gpartials;
  }

  // position of sun/moon, [m] in celestial RF
  Eigen::Matrix<double, 3, 1> rsun,rmon; 
  { // third body perturbations, Sun and Moon [m/sec^2] in celestial RF
    Eigen::Matrix<double, 3, 1> sun_acc;
    Eigen::Matrix<double, 3, 1> mon_acc;
    Eigen::Matrix<double, 3, 3> partials;
    dso::SunMoon(cmjd, r, params.GMSun, params.GMMon, sun_acc, mon_acc, rsun,
                 rmon, partials);
    f += (sun_acc + mon_acc);
    df += partials;
  }

  if (params.setide) { // earth tides on geopotential, gravity
    Eigen::Matrix<double, 3, 1> tacc;
    Eigen::Matrix<double, 3, 3> taccgrad;
    // Sun and Moon position in ECEF
    const Eigen::Matrix<double, 3, 1> rm_ecef = Rot.gcrf2itrf(rmon);
    const Eigen::Matrix<double, 3, 1> rs_ecef = Rot.gcrf2itrf(rsun);
      params.setide->acceleration(cmjd.tai2tt(), Rot.ut1(), r_geo, rm_ecef,
                                  rs_ecef, tacc, taccgrad);
    f += Rot.itrf2gcrf(tacc);
    const auto T = Rot.itrf2gcrf();
    df += T * taccgrad * T.transpose();
  }

  if (params.octide) 
  { // oean tides on geopotential, gravity
    Eigen::Matrix<double, 3, 1> tacc;
    params.octide->acceleration(cmjd.tai2tt(), r_geo, tacc);
    f += Rot.itrf2gcrf(tacc);
  }

  // Differential equation for the state transition matrix Φ(t, t_0) is:
  //            d/dt[Φ(t,t_0)] = F * Φ(t, t_0)
  // with initial conditions: Φ(t_0, t_0) = I_(6x6)
  // where F is the derivative of state transition matrix, aka
  //     |   0 (3x3)     I (3x3)   |
  // F = |                         | of size (6x6)
  //     | da/dr (3x3) da/dv (3x3) | 
  Eigen::Matrix<double, 6, 6> F; // dfdy
  {
    F.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Zero();
    F.block<3, 3>(0, 3) = Eigen::Matrix<double, 3, 3>::Identity();
    F.block<3, 3>(3, 0) = df;
    F.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Zero();
  }

  
  { // Derivative of combined state vector and state transition matrix
    // ---------------------------------------------------------------
    //             |   0 (3x3)     I (3x3)   |   | dr/dr0 (3x3) dr/dv0 (3x3) |
    // F*Φ(t,t0) = |                         | * |                           |
    //             | da/dr (3x3) da/dv (3x3) |   | dv/dr0 (3x3) dv/dv0 (3x3) |
    //       
    // |             dv/dr0                             dv/dv0               |
    // | (da/dr * dr/dr0 + da/dv * dv/dr0) (da/dr * dr/dv0 + da/dv * dv/dv0) |
    yPt.block<3, 1>(0, 0) = v;
    yPt.block<3, 1>(3, 0) = f;
    yPt.block<18, 1>(6, 0) = yP0.block<18, 1>(24, 0);

    Eigen::Matrix<double, 6, 6> Phi0;
    Phi0.block<1, 3>(0, 0) = yP0.block<3, 1>(6, 0).transpose();
    Phi0.block<1, 3>(1, 0) = yP0.block<3, 1>(9, 0).transpose();
    Phi0.block<1, 3>(2, 0) = yP0.block<3, 1>(12, 0).transpose();
    Phi0.block<1, 3>(0, 3) = yP0.block<3, 1>(15, 0).transpose();
    Phi0.block<1, 3>(1, 3) = yP0.block<3, 1>(18, 0).transpose();
    Phi0.block<1, 3>(2, 3) = yP0.block<3, 1>(21, 0).transpose();
    Phi0.block<1, 3>(3, 0) = yP0.block<3, 1>(24, 0).transpose();
    Phi0.block<1, 3>(4, 0) = yP0.block<3, 1>(27, 0).transpose();
    Phi0.block<1, 3>(5, 0) = yP0.block<3, 1>(30, 0).transpose();
    Phi0.block<1, 3>(3, 3) = yP0.block<3, 1>(33, 0).transpose();
    Phi0.block<1, 3>(4, 3) = yP0.block<3, 1>(36, 0).transpose();
    Phi0.block<1, 3>(5, 3) = yP0.block<3, 1>(39, 0).transpose();

    Eigen::Matrix<double, 3, 3> t1 =
        (F.block<3, 3>(3, 0) * Phi0.block<3, 3>(0, 0) +
         F.block<3, 3>(3, 3) * Phi0.block<3, 3>(3, 0))
            .transpose();
    yPt.block<9, 1>(24, 0) = Eigen::Map<Eigen::Matrix<double, 9, 1>>(
        t1.data(), t1.cols() * t1.rows());
    Eigen::Matrix<double, 3, 3> t2 =
        (F.block<3, 3>(3, 0) * Phi0.block<3, 3>(0, 3) +
         F.block<3, 3>(3, 3) * Phi0.block<3, 3>(3, 3))
            .transpose();
    yPt.block<9, 1>(33, 0) = Eigen::Map<Eigen::Matrix<double, 9, 1>>(
        t2.data(), t2.cols() * t2.rows());
  }

  return;
}

namespace {
//
//  yPhi0 = [ r, v, φ1, φ2, φ3 ]
//  where:
//  r  = [x y z]     : 3x1
//  v  = [Vx Vy Vz]  : 3x1
//  φ1 = dr/dy0      : 3xn
//       +----- 3 cols -------+-------- 3 cols -------+---- m cols ------|
//       | dx/x0 dx/dy0 dx/dz0 dx/dVx0 dx/dVy0 dx/dVz0 dx/dp1 dx/dp2 ... | -> index: [6,6+n)
//     = | dy/x0 dy/dy0 dy/dz0 dy/dVx0 dy/dVy0 dy/dVz0 dy/dp1 dy/dp2 ... | -> index: [6+n, 6+2n)
//       | dz/x0 dz/dy0 dz/dz0 dz/dVx0 dz/dVy0 dz/dVz0 dz/dp1 dz/dp2 ... | -> index: [6+2n, 6+3n)
// 
//  φ2 = dv/dy0      : 3xn
//  +----- 3 cols ----------+-------- 3 cols ----------+------ m cols ------|
//  | dVx/x0 dVx/dy0 dVx/dz0 dVx/dVx0 dVx/dVy0 dVx/dVz0 dVx/dp1 dVx/dp2 ... | -> index: [6+3n, 6+4n)
//  | dVy/x0 dVy/dy0 dVy/dz0 dVy/dVx0 dVy/dVy0 dVy/dVz0 dVy/dp1 dVy/dp2 ... | -> index: [6+4n, 6+5n)
//  | dVz/x0 dVz/dy0 dVz/dz0 dVz/dVx0 dVz/dVy0 dVz/dVz0 dVz/dp1 dVz/dp2 ... | -> index: [6+5n, 6+6n)
// 
//  φ3 = dp/dy0      : mxn
//  +----- 3 cols ----------+-------- 3 cols ----------+------ m cols ------|
//  | dp1/x0 dp1/dy0 dp1/dz0 dp1/dVx0 dp1/dVy0 dp1/dVz0 dp1/dp1 dp1/dp2 ... | -> index: [6+6n, 6+7n)
//  | dp2/x0 dp2/dy0 dp2/dz0 dp2/dVx0 dp2/dVy0 dp2/dVz0 dp2/dp1 dp2/dp2 ... | -> index: [6+7n, 6+8n)
//                                                                            -> index: [ ..., 6+(6+m)n )
//     = [ 0(mx6) I(mxm) ]
// 
Eigen::Matrix<double, m, n> phi3() noexcept {
  Eigen::Matrix<double, m, n> f3;
  f3.block<m, 6>(0, 0) = Eigen::Matrix<double, m, 6>::Zero();
  f3.block<m, m>(0, 6) = Eigen::Matrix<double, m, m>::Identity();
  return f3;
}
Eigen::Matrix<double, 3, n>
phi1(const Eigen::Matrix<double, (6 + 6*n), 1> &yp0) noexcept {
  Eigen::Matrix<double, 3, n> f1;
  constexpr const int of = 6;
  f1.block<1, n>(0, 0) = yp0.block<n, 1>(of, 0).transpose();
  f1.block<1, n>(1, 0) = yp0.block<n, 1>(of + n, 0).transpose();
  f1.block<1, n>(2, 0) = yp0.block<n, 1>(of + 2 * n, 0).transpose();
  return f1;
}
Eigen::Matrix<double, 3, n>
phi2(const Eigen::Matrix<double, (6 + 6*n), 1> &yp0) noexcept {
  Eigen::Matrix<double, 3, n> f2;
  constexpr const int of = 6 + 3 * n;
  f2.block<1, n>(0, 0) = yp0.block<n, 1>(of, 0).transpose();
  f2.block<1, n>(1, 0) = yp0.block<n, 1>(of + n, 0).transpose();
  f2.block<1, n>(2, 0) = yp0.block<n, 1>(of + 2 * n, 0).transpose();
  return f2;
}
} //unnamed namespace

void dso::VariationalEquations2(
    // seconds from reference epoch (TAI)
    double tsec,
    // state and state transition matrix (inertial RF) column-vector of size 6+6*n
    const Eigen::VectorXd &yP0,
    // state derivative and state transition matrix derivative (inertial RF)
    Eigen::Ref<Eigen::VectorXd> yPt, // column-vector of size 6+6*n
    // auxiliary parametrs
    dso::IntegrationParameters &params) noexcept {

  // current mjd, TAI
  const dso::TwoPartDate _cmjd(params.mjd_tai +
                              dso::TwoPartDate(0e0, tsec / 86400e0));
  const dso::TwoPartDate cmjd = _cmjd.normalized();

  // terrestrial to celestial for requested epoch t
  dso::Itrs2Gcrs Rot(cmjd.tai2tt(), &params.eopLUT);

  // split position and velocity vectors (GCRF) at t=t0
  Eigen::Matrix<double, 3, 1> r = yP0.block<3, 1>(0, 0);
  Eigen::Matrix<double, 3, 1> v = yP0.block<3, 1>(3, 0);

  // f = y''(t_0,y_0) = a(t_0, y_0) [=(a_x, a_y, a_z)]
  Eigen::Matrix<double, 3, 1> f  = Eigen::Matrix<double, 3, 1>::Zero();
  // df = grav(f(t_0,y_0)) = grad(a(t_0, y_0))  
  //                 | da_x/dx0 da_y/dy0 da_z/dz0 |
  //               = | da_x/dx0 da_y/dy0 da_z/dz0 |
  //                 | da_x/dx0 da_y/dy0 da_z/dz0 |
  Eigen::Matrix<double, 3, 3> dadr = Eigen::Matrix<double, 3, 3>::Zero();
  Eigen::Matrix<double, 3, 3> dadv = Eigen::Matrix<double, 3, 3>::Zero();
  Eigen::Matrix<double, m, 3> dadp = Eigen::Matrix<double, m, 3>::Zero();

  // satellite position at t=t0 (ITRF)
  Eigen::Matrix<double, 3, 1> r_itrf = Rot.gcrf2itrf(r);

  { // compute gravity-induced acceleration
    Eigen::Matrix<double, 3, 3> gradient;
    Eigen::Matrix<double, 3, 1> acc;
    dso::gravity_acceleration(params.harmonics, r_itrf, params.degree,
                   params.harmonics.Re(), params.harmonics.GM(), acc,
                   gradient, params.V, params.W);

    // transform acceleration and gradient to GCRF
    acc = Rot.itrf2gcrf(acc);
    const auto T = Rot.itrf2gcrf();
    gradient = T * gradient* T.transpose();

    f += acc;
    dadr += gradient;
  }

  // position of sun/moon, [m] in GCRF
  Eigen::Matrix<double, 3, 1> rsun,rmon; 
  { // third body perturbations, Sun and Moon [m/sec^2] in celestial RF
    Eigen::Matrix<double, 3, 1> sun_acc;
    Eigen::Matrix<double, 3, 1> mon_acc;
    Eigen::Matrix<double, 3, 3> gradient;
    dso::SunMoon(cmjd, r, params.GMSun, params.GMMon, sun_acc, mon_acc, rsun,
                 rmon, gradient);
    f += (sun_acc + mon_acc);
    dadr += gradient;
  }

  if (params.setide) { 
    // earth tides on geopotential
    Eigen::Matrix<double, 3, 1> acc;
    Eigen::Matrix<double, 3, 3> gradient;
    // Sun and Moon position in ECEF
    const Eigen::Matrix<double, 3, 1> rm_ecef = Rot.gcrf2itrf(rmon);
    const Eigen::Matrix<double, 3, 1> rs_ecef = Rot.gcrf2itrf(rsun);
      params.setide->acceleration(cmjd.tai2tt(), Rot.ut1(), r_itrf, rm_ecef,
                                  rs_ecef, acc, gradient);
    // acceleration and gradient to GCRF
    f += Rot.itrf2gcrf(acc);
    const auto T = Rot.itrf2gcrf();
    dadr += T * gradient * T.transpose();
  }

  if (params.octide) {
    // oean tides on geopotential, gravity
    Eigen::Matrix<double, 3, 1> acc;
    params.octide->acceleration(cmjd.tai2tt(), r_itrf, acc);
    f += Rot.itrf2gcrf(acc);
  }
  
  Eigen::Matrix<double,n,n> f123;
  {
    f123.block<3,n>(0,0) = phi1(yP0);
    f123.block<3,n>(3,0) = phi2(yP0);
    if (m) f123.block<m,n>(6,0) = phi3();
  }

  Eigen::Matrix<double,n,n> A;
  {
    A.block<3,3>(0,0) = Eigen::Matrix<double, 3, 3>::Zero();
    A.block<3,3>(0,3) = Eigen::Matrix<double, 3, 3>::Identity();
    if (m) A.block<3,m>(0,6) = Eigen::Matrix<double, 3, m>::Zero();

    A.block<3,3>(3,0) = dadr;
    A.block<3,3>(3,3) = dadv;
    if (m) A.block<3,m>(3,6) = dadp;

    if (m) {
    A.block<m,3>(6,0) = Eigen::Matrix<double, m, 3>::Zero();
    A.block<m,3>(6,3) = Eigen::Matrix<double, m, 3>::Zero();
    A.block<m,m>(6,6) = Eigen::Matrix<double, m, m>::Zero();
    }
  }
  Eigen::Matrix<double,n,n> F = A * f123;

  // Assign to yPt
  yPt.block<3,1>(0,0) = v;
  yPt.block<3,1>(3,0) = f;
  yPt.block<n,1>(6+0*n,0) = F.block<1,n>(0,0).transpose(); // dφ1(t,t0) / dt
  yPt.block<n,1>(6+1*n,0) = F.block<1,n>(1,0).transpose();
  yPt.block<n,1>(6+2*n,0) = F.block<1,n>(2,0).transpose();
  yPt.block<n,1>(6+3*n,0) = F.block<1,n>(3,0).transpose(); // dφ2(t,t0) / dt
  yPt.block<n,1>(6+4*n,0) = F.block<1,n>(4,0).transpose();
  yPt.block<n,1>(6+5*n,0) = F.block<1,n>(5,0).transpose();

  return;
}
