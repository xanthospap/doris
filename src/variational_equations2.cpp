#include "astrodynamics.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "orbit_integration.hpp"
#include "iers2010/cel2ter.hpp"
#include <datetime/dtcalendar.hpp>
#include <geodesy/ellipsoid.hpp>

constexpr const int m = 0;
constexpr const int F1 = 3 * (6+m);
constexpr const int F2 = 3 * (6+m);
constexpr const int F3 = m * (6+m);

void dso::VariationalEquations(
    /* seconds from reference epoch (TAI) */
    double tsec,
    /* state and state transition matrix (inertial RF) column-vector */
    const Eigen::VectorXd &yP0,
    /* state derivative and state transition matrix derivative (inertial RF) */
    Eigen::Ref<Eigen::VectorXd> yPt,
    /* auxiliary parametrs */
    dso::IntegrationParameters *params) noexcept
{
  /*          yP0 =                         yPt
   *    --------------------------------------------------------------------
   *     r     (3x1)   -> 3                 v
   *     v     (3x1)   -> 3                 a
   *     φ1    (3x6+m) -> 3 * (6+m)         φ1'
   *     φ2    (3x6+m) -> 3 * (6+m)         φ2'
   *     φ3    (mx6+m) -> m * (6+m)         φ3' (==0)
   *     -------------------------------------------------------------------
   *                      6 + 2*(3 * (6+m)) + m * (6+m)
   */

  assert( yP0.rows() == 6 + 2*(3 * (6+m)) + m * (6+m) );
  assert( yP0.cols() == 1 );
  
  /* current mjd, TAI */
  const dso::TwoPartDate cmjd(params->mjd_tai +
                              dso::TwoPartDate(0e0, tsec / 86400e0));

  /* terrestrial to celestial for requested epoch t */
  dso::Itrs2Gcrs Rot(cmjd.tai2tt(), &(params->eop_lookup_table()));

  /* split position and velocity vectors (GCRF) at t=t0 */
  Eigen::Matrix<double, 3, 1> r = yP0.block<3, 1>(0, 0);
  Eigen::Matrix<double, 3, 1> v = yP0.block<3, 1>(3, 0);

  Eigen::Matrix<double, 3, 1> a    = Eigen::Matrix<double, 3, 1>::Zero();
  Eigen::Matrix<double, 3, 3> dadr = Eigen::Matrix<double, 3, 3>::Zero();
  Eigen::Matrix<double, 3, 3> dadv = Eigen::Matrix<double, 3, 3>::Zero();

  /* satellite position at t=t0 (ITRF) */
  Eigen::Matrix<double, 3, 1> r_itrf = Rot.gcrf2itrf(r);

  { /* compute gravity-induced acceleration */
    Eigen::Matrix<double, 3, 3> gradient;
    Eigen::Matrix<double, 3, 1> acc;
    if (params->egravity->acceleration(r_itrf, acc, gradient)) {
      fprintf(stderr,
              "[ERROR] Failed computing gravity/gradient (traceback: %s)\n",
              __func__);
    }

    /* transform acceleration and gradient to GCRF */
    acc = Rot.itrf2gcrf(acc);
    const auto T = Rot.itrf2gcrf();
    gradient = T * gradient* T.transpose();

    a += acc;
    dadr += gradient;
  }

  /* yPt(0:3) = v */
  yPt.block<3,1>(0,0) = v;
  /* yPt(3:6) = a */
  yPt.block<3,1>(3,0) = a;
  /* φ1'(t,t0) = φ2(t,t0) */
  yPt.block<F1,1>(6,0) = yP0.block<F1,1>(6+F1,0);
  /* φ2'(t,t0) = da/dr*φ1(t,t0) + da/dv*φ2(t,t0) + da/dp*φ3(t,t0) */
  yPt.block<F2, 1>(6 + F1, 0) =
      (dadr * yP0.block<F1, 1>(6, 0).reshaped(6 + m, 3).transpose() +
       dadv * yP0.block<F2, 1>(6 + F1, 0).reshaped(6 + m, 3).transpose())
          .transpose()
          .reshaped();
  /* φ3'(t,t0) = 0 */

  return;
}
