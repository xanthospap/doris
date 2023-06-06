#include "astrodynamics.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "iers2010/cel2ter.hpp"
#include "orbit_integration.hpp"

constexpr const int m = 0;
constexpr const int F1 = 3 * (6 + m);
constexpr const int F2 = 3 * (6 + m);
constexpr const int F3 = m * (6 + m);

void dso::VariationalEquations(
    /* seconds from reference epoch (TAI) */
    double tsec,
    /* state and state transition matrix (inertial RF) column-vector */
    const Eigen::VectorXd &yP0,
    /* state derivative and state transition matrix derivative (inertial RF) */
    Eigen::Ref<Eigen::VectorXd> yPt,
    /* auxiliary parametrs */
    dso::IntegrationParameters *params) noexcept {
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

  assert(yP0.rows() == 6 + 2 * (3 * (6 + m)) + m * (6 + m));
  assert(yP0.cols() == 1);
  auto error = dso::iStatus(0);

  /* current mjd, TAI */
  auto cmjd = params->tai();
  cmjd.add_seconds(tsec);

  /* terrestrial to celestial for requested epoch t */
  dso::Itrs2Gcrs Rot(cmjd.tai2tt(), &(params->eop_lookup_table()));

  /* split position and velocity vectors (GCRF) at t=t0 */
  Eigen::Matrix<double, 3, 1> r = yP0.block<3, 1>(0, 0);
  Eigen::Matrix<double, 3, 1> v = yP0.block<3, 1>(3, 0);

  Eigen::Matrix<double, 3, 1> a = Eigen::Matrix<double, 3, 1>::Zero();
  Eigen::Matrix<double, 3, 3> dadr = Eigen::Matrix<double, 3, 3>::Zero();
  Eigen::Matrix<double, 3, 3> dadv = Eigen::Matrix<double, 3, 3>::Zero();

  /* satellite position at t=t0 (ITRF) */
  Eigen::Matrix<double, 3, 1> r_itrf = Rot.gcrf2itrf(r);

  /* Sun and Moon position in GCRF */
  Eigen::Matrix<double, 3, 1> sun_gcrf, sun_itrf;
  Eigen::Matrix<double, 3, 1> mon_gcrf, mon_itrf;
  {
    error += dso::planet_pos(dso::Planet::SUN, Rot.tt(), sun_gcrf);
    error += dso::planet_pos(dso::Planet::MOON, Rot.tt(), mon_gcrf);
    sun_itrf = Rot.gcrf2itrf(sun_gcrf);
    mon_itrf = Rot.gcrf2itrf(mon_gcrf);
    if (error) {
      fprintf(stderr,
              "[ERROR] Failed computing planetary ephemeris (traceback: %s)\n",
              __func__);
    }
  }

  {
    /* compute geopotential coefficients for solid earth tides */
    params->solid_earth_tide()->operator()(Rot.tt(), Rot.ut1(), mon_itrf, sun_itrf);
  }

  {
    /* compute geopotential coefficients for ocean tides */
    if (params->ocean_tide()->operator()(Rot.tt(), Rot.ut1())) {
      error += dso::iStatus(1);
      fprintf(stderr,
              "[ERROR] Failed computing ocean tides geopotential corrections "
              "(traceback: %s)\n",
              __func__);
    }
  }

  {
    /* compute geopotential coefficients for ocean pole tides */
    if (params->ocean_pole_tide()->operator()(Rot.tt(), Rot.eop().xp, Rot.eop().yp)) {
      error += dso::iStatus(1);
      fprintf(stderr,
              "[ERROR] Failed computing ocean tides geopotential corrections "
              "(traceback: %s)\n",
              __func__);
    }
    /* solid earth pole tide */
    [[maybe_unused]]const auto dcs =
        params->solid_earth_pole_tide()->delta_stokes_21(Rot.tt(), Rot.eop().xp, Rot.eop().yp);
  }

  /* geopotential coefficients (gravity+tides) to acceleration */
  {
    Eigen::Matrix<double, 3, 3> gradient;
    Eigen::Matrix<double, 3, 1> acc;
    dso::gravity_acceleration(
        params->accumulate_geopotential_coeffs(), r_itrf,
        params->earth_gravity()->max_degree(),
        params->earth_gravity()->geopotential_coeffs().Re(),
        params->earth_gravity()->geopotential_coeffs().GM(), acc, gradient,
        &(params->workspace_v()), &(params->workspace_w()));
    /* transform acceleration and gradient to GCRF */
    acc = Rot.itrf2gcrf(acc);
    const auto T = Rot.itrf2gcrf();
    gradient = T * gradient * T.transpose();
    /* add to global matrices */
    a += acc;
    dadr += gradient;
  }

  /* acceleration (correction) due to relativity, IERS-2010 */
  {
    /* Time derivative of the position of the Earth with respect to Sun */
    auto dRes = sun_gcrf;
    {
      const double dt = 10e0; /* seconds */
      auto t2 = cmjd.tai2tt();
      t2.add_seconds(dt);
      Eigen::Matrix<double, 3, 1> sun2;
      if (dso::planet_pos(dso::Planet::SUN, t2, sun2)) {
        fprintf(stderr, "ERROR Failed extracting Sun/Moon position!\n");
        return;
      }
      dRes = (sun2 - sun_gcrf) / dt;
    }
    a += iers2010::relativistic_correction(r, sun_gcrf, dRes, v);
  }

  printf("Acceleration a = %.9e %.9e %.9e\n", a(0), a(1), a(2));

  /* (direct) solar radiation pressure */
  {
    /* shadow factor */
    const double f = dso::conic_shadow_factor(r, sun_gcrf, 6.957e8);
    if (f!=0e0) {
      /* get macromodel in ECI */
      std::vector<dso::MacroModelComponent> sv;
      if (params->svframe()->macromodel_iterator(cmjd,sv)) {
        fprintf(stderr, "[ERROR] Failed getting attitude (traceback: %s)\n",
                __func__);
      }
      /* compute SRP */
      const Eigen::Matrix<double, 3, 1> asrp =
          dso::direct_solar_radiation_pressure(sv, r, sun_gcrf,
                                               params->svframe()->mass());
      auto Asrp = asrp * f;
      printf("SRP a = %.9e %.9e %.9e\n", Asrp(0), Asrp(1), Asrp(2));
    }
  }

  /* yPt(0:3) = v */
  yPt.block<3, 1>(0, 0) = v;
  /* yPt(3:6) = a */
  yPt.block<3, 1>(3, 0) = a;
  /* φ1'(t,t0) = φ2(t,t0) */
  yPt.block<F1, 1>(6, 0) = yP0.block<F1, 1>(6 + F1, 0);
  /* φ2'(t,t0) = da/dr*φ1(t,t0) + da/dv*φ2(t,t0) + da/dp*φ3(t,t0) */
  yPt.block<F2, 1>(6 + F1, 0) =
      (dadr * yP0.block<F1, 1>(6, 0).reshaped(6 + m, 3).transpose() +
       dadv * yP0.block<F2, 1>(6 + F1, 0).reshaped(6 + m, 3).transpose())
          .transpose()
          .reshaped();
  /* φ3'(t,t0) = 0 */

  return;
}
