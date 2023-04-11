#include "astrodynamics.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "orbit_integration.hpp"
#include "iers2010/cel2ter.hpp"
#include <datetime/dtcalendar.hpp>
#include <exception>
#include <thread>

constexpr const int accountforpoletide = true;
constexpr const int m = 2;
constexpr const int n = 6 + m;

namespace {
void gravity(const Eigen::Matrix<double, 3, 1> &ritrf,
             const dso::Itrs2Gcrs &Rot,
             const dso::IntegrationParameters &params,
             Eigen::Matrix<double, 3, 1> &acc,
             Eigen::Matrix<double, 3, 3> &dadr,
             Eigen::Matrix<double, 3, 3> &dadv, int &error) noexcept {
  try {
  /* preset error */
  error = 0;

  /* acceleration and gradient in ITRF */
  Eigen::Matrix<double, 3, 1> a;
  Eigen::Matrix<double, 3, 3> gradient;
  test::gravacc3(params.harmonics, ritrf, params.degree, params.harmonics.Re(),
                 params.harmonics.GM(), a, gradient, params.V, params.W);

  /* transform acceleration and gradient to GCRF */
  acc = Rot.itrf2gcrf(a);
  const auto T = Rot.itrf2gcrf();
  dadr = T * gradient * T.transpose();

  /* gradient wrt velocity is zero */
  dadv = Eigen::Matrix<double, 3, 3>::Zero();
  } catch (std::exception &e) {
    fprintf(stderr, "[ERROR] exception caught in %s; type is %s\n", __func__, e.what());
    error = 1;
  }

  return;
}

[[maybe_unused]]
double sv_in_shadow(const Eigen::Matrix<double, 3, 1> &rgcrf,
                 const Eigen::Matrix<double, 3, 1> &sun_gcrf,
                 double Re, double Rs=6.957e8) noexcept {
  // return rgcrf.dot(sun_gcrf) < -std::sqrt(rgcrf.squaredNorm() - Re * Re);

  const double th_e = std::asin(Re / rgcrf.norm()); // b
  const double th_s = std::asin(Rs / (sun_gcrf - rgcrf).norm()); // a
  const double th_es = std::acos(-rgcrf.dot(sun_gcrf - rgcrf) /
                                 (rgcrf.norm() * (sun_gcrf - rgcrf).norm()));

  /* check if we are in the penumbra area */
  if ((std::abs(th_s-th_e) < th_es) && (th_es < th_s+th_e)) {
    fprintf(stderr, "[LOG] SV in partial shadow\n");

    const double caf = std::acos((th_s * th_s + th_es * th_es - th_e * th_e) /
                                 (2e0 * th_s * th_es));
    const double cbd = std::acos((th_e * th_e + th_es * th_es - th_s * th_s) /
                                 (2e0 * th_e * th_es));
    const double S_afc = 0.5e0 * caf * th_s * th_s;
    const double S_aec = 0.5e0 * (th_s * std::sin(caf)) * (th_s * std::cos(caf));
    const double S_bdc = 0.5e0 * cbd * th_e * th_e;
    const double S_bec = 0.5e0 * (th_e * std::sin(cbd)) * (th_e * std::cos(cbd));

    const double S = 2e0 * (S_afc - S_aec) + 2e0 * (S_bdc - S_bec);
    return std::min(1e0 - S / iers2010::DPI / (th_s*th_s), 1e0);
  } else {
    return (th_s+th_e <= th_es); /* no occultation */
  }
}

[[maybe_unused]]
double sv_in_shadow_conic(const Eigen::Matrix<double, 3, 1> &rgcrf,
                 const Eigen::Matrix<double, 3, 1> &sun_gcrf,
                 double Re) noexcept {
  return rgcrf.dot(sun_gcrf) < -std::sqrt(rgcrf.squaredNorm() - Re * Re);
}

int solar_radiation_pressure(const dso::TwoPartDate &tai,
                             const Eigen::Matrix<double, 3, 1> &rgcrf,
                             const Eigen::Matrix<double, 3, 1> &sun_gcrf,
                             const Eigen::Quaternion<double> q,
                             dso::IntegrationParameters &params, double Cr,
                             Eigen::Matrix<double, 3, 1> &acc,
                             Eigen::Matrix<double, 3, 3> &dadr,
                             Eigen::Matrix<double, 3, 3> &dadv,
                             Eigen::Matrix<double, 3, 1> &dadC) {
  static double oc_factor = 0;
  double ShadowFactor = sv_in_shadow(rgcrf, sun_gcrf, iers2010::Re);
  if (oc_factor != ShadowFactor) {
    fprintf(stderr, "[WRNNG] Changing shadow status from %.3f to %.3f\n",
            oc_factor, ShadowFactor);
    printf("[SHDW ] Changing shadow status from %.3f to %.3f at %.9f\n",
            oc_factor, ShadowFactor, tai.mjd());
    oc_factor = ShadowFactor;
  }

  // sun solar radiation pressure, Montenbruck, Sec 3.4
  constexpr const double Ps = 4.56e-6; // Nm^(-2)

  // unit vector directed from spacecrafe to Sun (GCRF)
  const Eigen::Matrix<double, 3, 1> svsun = (sun_gcrf - rgcrf).normalized();

  // get attitude rotation matrix / quaternion
  //Eigen::Quaternion<double> q;
  //assert(!params.svFrame->get_attitude_quaternion(tai, q));

  // loop over box-wings
  const int numPlates = params.svFrame->plates().NumPlates;
  const dso::MacroModelComponent *plate = nullptr;
  Eigen::Matrix<double, 3, 1> S = Eigen::Matrix<double, 3, 1>::Zero();
  for (int i = 0; i < numPlates; i++) {
    plate = params.svFrame->plates().mmcomponents + i;
    /* (unit) vector of plate, normal in sv-fixed RF */
    Eigen::Matrix<double, 3, 1> n(plate->m_normal);
    /* inclination of the ith plate to the SV-to-Sun vector */
    const auto ni = q.conjugate().normalized() * n;
    const double cosA = ni.transpose() * svsun;
    if (cosA > 0e0) {
      const double *R = plate->m_optical_properties;
      S += (plate->m_surf * cosA) *
           (2e0 * (R[1] / 3e0 + R[0] * cosA) * ni + (1e0 - R[0]) * svsun);
    }
  }
  S /= params.svFrame->mass();

  // scale
  const double rss = (sun_gcrf - rgcrf).norm();
  const double fac = -Cr * ShadowFactor * Ps *
                     std::pow(iers2010::AU / rss, 2e0);

  // acceleration
  acc = fac * S;
  // partials
  const auto rs = rgcrf - sun_gcrf;
  dadr = -fac *
         (Eigen::Matrix<double, 3, 3>::Identity() -
          3e0 * rs * rs.transpose() / rss / rss) /
         rss;
  dadv = Eigen::Matrix<double, 3, 3>::Zero();
  dadC = acc / Cr;

  return 0;
}

void tides_n_tba(const dso::TwoPartDate &tai,
                 const Eigen::Matrix<double, 3, 1> &ritrf,
                 const Eigen::Matrix<double, 3, 1> &rgcrf,
                 const Eigen::Matrix<double, 3, 1> &sungcrf,
                 const Eigen::Matrix<double, 3, 1> &moongcrf,
                 const dso::Itrs2Gcrs &Rot,
                 const dso::IntegrationParameters &params,
                 Eigen::Matrix<double, 3, 1> &acc,
                 Eigen::Matrix<double, 3, 3> &dadr,
                 Eigen::Matrix<double, 3, 3> &dadv, int &error) noexcept {
  error = 0;

  /* third body attraction, point mass */
  {
    Eigen::Matrix<double, 3, 3> sgrad, mgrad;
    const double GMSun = params.GMSun;
    const double GMMoon = params.GMMon;
    const Eigen::Matrix<double, 3, 1> sun_acc =
        dso::point_mass_accel(GMSun * 1e9, rgcrf, sungcrf, sgrad);
    const Eigen::Matrix<double, 3, 1> moon_acc =
        dso::point_mass_accel(GMMoon * 1e9, rgcrf, moongcrf, mgrad);
    acc = (sun_acc + moon_acc);
    dadr = (sgrad + mgrad);
  }

  /* solid earth tides */
  if (params.setide) {
    Eigen::Matrix<double, 3, 1> tacc;
    Eigen::Matrix<double, 3, 3> tgradient;
    /* Sun and Moon position in ITRF */
    const Eigen::Matrix<double, 3, 1> rm_ecef = Rot.gcrf2itrf(moongcrf);
    const Eigen::Matrix<double, 3, 1> rs_ecef = Rot.gcrf2itrf(sungcrf);
    if (accountforpoletide) {
      params.setide->acceleration(tai.tai2tt(), Rot.ut1(), Rot.eop().xp,
                                  Rot.eop().yp, ritrf, rm_ecef, rs_ecef, tacc,
                                  tgradient);
    } else {
      Eigen::Matrix<double, 3, 1> pacc;
      Eigen::Matrix<double, 3, 3> pgradient;
      params.setide->acceleration(tai.tai2tt(), Rot.ut1(), ritrf, rm_ecef,
                                  rs_ecef, tacc, tgradient);
    }
    // acceleration and gradient to GCRF
    acc += Rot.itrf2gcrf(tacc);
    const auto T = Rot.itrf2gcrf();
    dadr += T * tgradient * T.transpose();
  }

  /* ocean tides */
  if (params.octide) {
    Eigen::Matrix<double, 3, 1> tacc;
    params.octide->acceleration(tai.tai2tt(), ritrf, tacc);
    acc += Rot.itrf2gcrf(tacc);
  }

  /* gradient wrt satellite velocity */
  dadv = Eigen::Matrix<double, 3, 3>::Zero();

  return;
}

void drag(const dso::TwoPartDate &tai, const Eigen::Matrix<double, 3, 1> &rgcrf,
          const Eigen::Matrix<double, 3, 1> &vgcrf,
          const Eigen::Matrix<double, 3, 1> &ritrf, const dso::Itrs2Gcrs &Rot,
          const Eigen::Quaternion<double> &q,
          dso::IntegrationParameters &params, Eigen::Matrix<double, 3, 1> &acc,
          Eigen::Matrix<double, 3, 3> &dadr, Eigen::Matrix<double, 3, 3> &dadv,
          Eigen::Matrix<double, 3, 1> &dadC, int &error) noexcept {
  /* set error to noerror */
  error = 0;
  
  try {
  /* get attitude rotation matrix / quaternion */
  // Eigen::Quaternion<double> q;
  //if (params.svFrame->get_attitude_quaternion(tai, q)) {
  //  error = 1;
  //  return;
  //}

  /* relative velocity wrt atmosphere */
  Eigen::Matrix<double, 3, 1> vrel;
  {
    // Velocity relative to the Earth's atmosphere
    const Eigen::Matrix<double, 3, 1> omega{0e0, 0e0, Rot.omega_earth()};
    vrel = vgcrf - omega.cross(rgcrf);
  }

  /* loop over box-wings to get projected area */
  double S = 0e0;
  const int numPlates = params.svFrame->plates().NumPlates;
  const dso::MacroModelComponent *plate = nullptr;
  for (int i = 0; i < numPlates; i++) {
    plate = params.svFrame->plates().mmcomponents + i;
    /* (unit) vector of plate, normal in sv-fixed RF */
    Eigen::Matrix<double, 3, 1> n(plate->m_normal);
    /* inclination of the ith plate to the  vector */
    const auto ni = q.conjugate().normalized() * n;
    const double cosA = ni.transpose() * vrel.normalized();
    if (cosA > 0e0) {
      S += (plate->m_surf * cosA);
    }
  }

  /* prjected area / satellite mass */
  S /= params.svFrame->mass();

  /* get atmospheric density, using the UTC date */
  const dso::TwoPartDate utc = tai.tai2utc();
  const auto lfh = dso::car2ell<dso::ellipsoid::grs80>(ritrf);
  if (params.Dtm20.set(utc, lfh[0], lfh[1], lfh[2] / 1e3)) {
    error = 2;
    return;
  }
  if (params.Dtm20.dtm3()) {
    error = 2;
    return;
  }
  const double rho = params.Dtm20.total_density_kgm3(); /* density here */

  /* get dρ/dr aka gradient of density w.r.t position */
  Eigen::Matrix<double, 3, 1> drhodr;
  {
    const double dr = 5e0;
    // w.r.t X component
    auto dlfh = dso::car2ell<dso::ellipsoid::grs80>(
        ritrf + dr * Eigen::Matrix<double, 3, 1>({1e0, 0e0, 0e0}));
    error += params.Dtm20.set(utc, dlfh[0], dlfh[1], dlfh[2] / 1e3);
    error += params.Dtm20.dtm3();
    double rhop = params.Dtm20.total_density_grcm3();
    dlfh = dso::car2ell<dso::ellipsoid::grs80>(
        ritrf - dr * Eigen::Matrix<double, 3, 1>({1e0, 0e0, 0e0}));
    error += params.Dtm20.set(utc, dlfh[0], dlfh[1], dlfh[2] / 1e3);
    error += params.Dtm20.dtm3();
    double rhom = params.Dtm20.total_density_grcm3();
    drhodr(0) = (rhop - rhom) / (2 * dr);
    // w.r.t Y component
    dlfh = dso::car2ell<dso::ellipsoid::grs80>(
        ritrf + dr * Eigen::Matrix<double, 3, 1>({0e0, 1e0, 0e0}));
    error += params.Dtm20.set(utc, dlfh[0], dlfh[1], dlfh[2] / 1e3);
    error += params.Dtm20.dtm3();
    rhop = params.Dtm20.total_density_grcm3();
    dlfh = dso::car2ell<dso::ellipsoid::grs80>(
        ritrf - dr * Eigen::Matrix<double, 3, 1>({0e0, 1e0, 0e0}));
    error += params.Dtm20.set(utc, dlfh[0], dlfh[1], dlfh[2] / 1e3);
    error += params.Dtm20.dtm3();
    rhom = params.Dtm20.total_density_grcm3();
    drhodr(1) = (rhop - rhom) / (2 * dr);
    // w.r.t Z component
    dlfh = dso::car2ell<dso::ellipsoid::grs80>(
        ritrf + dr * Eigen::Matrix<double, 3, 1>({0e0, 0e0, 1e0}));
    error += params.Dtm20.set(utc, dlfh[0], dlfh[1], dlfh[2] / 1e3);
    error += params.Dtm20.dtm3();
    rhop = params.Dtm20.total_density_grcm3();
    dlfh = dso::car2ell<dso::ellipsoid::grs80>(
        ritrf - dr * Eigen::Matrix<double, 3, 1>({0e0, 0e0, 1e0}));
    error += params.Dtm20.set(utc, dlfh[0], dlfh[1], dlfh[2] / 1e3);
    error += params.Dtm20.dtm3();
    rhom = params.Dtm20.total_density_grcm3();
    drhodr(2) = (rhop - rhom) / (2 * dr);
    // scale (gram/cm3 to kg/m^3)
    drhodr *= 1e-7;
  }

  /* compute drag acceleration and gradients */
  acc = dso::drag_accel(vrel, S, params.Cd, rho, drhodr, dadr, dadv, dadC);
  } catch (std::exception &e) {
    fprintf(stderr, "[ERROR] exception caught in %s; type is %s\n", __func__, e.what());
    error = 2;
  }

  /* all done */
  return;
}

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
}//unnamed namespace

void dso::VariationalEquations_thread(
    // seconds from reference epoch (TAI)
    double tsec,
    // state and state transition matrix (inertial RF) column-vector of size 6+6*n
    const Eigen::VectorXd &yP0,
    // state derivative and state transition matrix derivative (inertial RF)
    Eigen::Ref<Eigen::VectorXd> yPt, // column-vector of size 6+6*n
    // auxiliary parametrs
    dso::IntegrationParameters &params) noexcept {

  /* current mjd, TAI */
  const dso::TwoPartDate _cmjd(params.mjd_tai +
                              dso::TwoPartDate(0e0, tsec / 86400e0));
  const dso::TwoPartDate cmjd = _cmjd.normalized();

  /* terrestrial to celestial for requested epoch t */
  dso::Itrs2Gcrs Rot(cmjd.tai2tt(), &params.eopLUT);

  /* split position and velocity vectors (GCRF) at t=t0 */
  Eigen::Matrix<double, 3, 1> r = yP0.block<3, 1>(0, 0);
  Eigen::Matrix<double, 3, 1> v = yP0.block<3, 1>(3, 0);

  Eigen::Matrix<double, 3, 1> a    = Eigen::Matrix<double, 3, 1>::Zero();
  Eigen::Matrix<double, 3, 3> dadr = Eigen::Matrix<double, 3, 3>::Zero();
  Eigen::Matrix<double, 3, 3> dadv = Eigen::Matrix<double, 3, 3>::Zero();
  Eigen::Matrix<double, 3, m> dadp = Eigen::Matrix<double, 3, m>::Zero();

  /* satellite position at t=t0 (ITRF) */
  Eigen::Matrix<double, 3, 1> r_itrf = Rot.gcrf2itrf(r);

  /* Get the attitude matrix so that we don't try to get it within the threads 
   * and mess with the stream
   */
  
  /* gravity acceleration and gradients */
  Eigen::Matrix<double, 3, 1> gravity_a;
  Eigen::Matrix<double, 3, 3> gravity_dadr;
  Eigen::Matrix<double, 3, 3> gravity_dadv;
  int gravity_error = 0;
  std::thread tgravity(gravity, std::cref(r_itrf), std::cref(Rot),
                       std::cref(params), std::ref(gravity_a),
                       std::ref(gravity_dadr), std::ref(gravity_dadv),
                       std::ref(gravity_error));
  
  /* position of Sun and Moon in GCRF at this instance */
  Eigen::Matrix<double, 3, 1> rSun,rMon;
  {
    const auto mjd_tt = cmjd.tai2tt();
    const double jd = (mjd_tt._big + dso::mjd0_jd) + mjd_tt._small;
    double rsun[3], rmon[3];

    /* position vector of sun/moon, in J2000, [km] */
    dso::cspice::j2planet_pos_from(dso::cspice::jd2et(jd), 10, 399, rsun);
    dso::cspice::j2planet_pos_from(dso::cspice::jd2et(jd), 301, 399, rmon);

    rSun = Eigen::Matrix<double, 3, 1>(rsun); /* [m] */
    rMon = Eigen::Matrix<double, 3, 1>(rmon); /* [m] */

    rSun *= 1e3; /* [m] */
    rMon *= 1e3; /* [m] */
  }

  /* tides and third body acceleration and gradients */
  Eigen::Matrix<double, 3, 1> tidestba_a;
  Eigen::Matrix<double, 3, 3> tidestba_dadr;
  Eigen::Matrix<double, 3, 3> tidestba_dadv;
  int tidestba_error = 0;
  std::thread ttidestba(tides_n_tba, std::cref(cmjd), std::cref(r_itrf),
                        std::cref(r), std::cref(rSun), std::cref(rMon),
                        std::cref(Rot), std::cref(params), std::ref(tidestba_a),
                        std::ref(tidestba_dadr), std::ref(tidestba_dadv),
                        std::ref(tidestba_error));

    /* drag acceleration and gradients*/
    Eigen::Matrix<double, 3, 1> drag_a;
    Eigen::Matrix<double, 3, 3> drag_dadr;
    Eigen::Matrix<double, 3, 3> drag_dadv;
    Eigen::Matrix<double, 3, 1> drag_dadC;
    int drag_error;
    std::thread tdrag(drag, cmjd.tai2tt(), std::cref(r), std::cref(v),
                      std::cref(r_itrf), std::cref(Rot), std::ref(q),
                      std::ref(params), std::ref(drag_a), std::ref(drag_dadr),
                      std::ref(drag_dadv), std::ref(drag_dadC),
                      std::ref(drag_error));

    /* solar radiation pressure in main thread */
    {
      Eigen::Matrix<double, 3, 1> ta;
      Eigen::Matrix<double, 3, 3> tdadr;
      Eigen::Matrix<double, 3, 3> tdadv;
      Eigen::Matrix<double, 3, 1> tdadp;
      assert(!solar_radiation_pressure(cmjd, r, rSun, q, params, params.Cr, ta,
                                       tdadr, tdadv, tdadp));
      a = ta;
      dadr = tdadr;
      dadv = tdadv;
      dadp.block<3, 1>(0, 1) = tdadp;
    }

    Eigen::Matrix<double, n, n> f123;
    {
      f123.block<3, n>(0, 0) = phi1(yP0);
      f123.block<3, n>(3, 0) = phi2(yP0);
      if (m)
        f123.block<m, n>(6, 0) = phi3();
    }

    /* join threads; from here on we need the computations */
    ttidestba.join();
    tgravity.join();
    tdrag.join();

    /* add acceleration and gradients */
    a += (gravity_a + drag_a + tidestba_a);
    dadr += (gravity_dadr + drag_dadr + tidestba_dadr);
    dadv += (gravity_dadv + drag_dadv + tidestba_dadv);
    dadp.block<3, 1>(0, 0) = drag_dadC;
    assert(!(gravity_error + drag_error + tidestba_error));

    /* construct system of equations */
    Eigen::Matrix<double, n, n> A;
    {
      A.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Zero();
      A.block<3, 3>(0, 3) = Eigen::Matrix<double, 3, 3>::Identity();
      if (m)
        A.block<3, m>(0, 6) = Eigen::Matrix<double, 3, m>::Zero();

      A.block<3, 3>(3, 0) = dadr;
      A.block<3, 3>(3, 3) = dadv;
      if (m)
        A.block<3, m>(3, 6) = dadp;

      if (m) {
        A.block<m, 3>(6, 0) = Eigen::Matrix<double, m, 3>::Zero();
        A.block<m, 3>(6, 3) = Eigen::Matrix<double, m, 3>::Zero();
        A.block<m, m>(6, 6) = Eigen::Matrix<double, m, m>::Zero();
      }
    }
    Eigen::Matrix<double, n, n> F = A * f123;

    // Assign to yPt
    yPt.block<3, 1>(0, 0) = v;
    yPt.block<3, 1>(3, 0) = a;
    yPt.block<n, 1>(6 + 0 * n, 0) =
        F.block<1, n>(0, 0).transpose(); // dφ1(t,t0) / dt
    yPt.block<n, 1>(6 + 1 * n, 0) = F.block<1, n>(1, 0).transpose();
    yPt.block<n, 1>(6 + 2 * n, 0) = F.block<1, n>(2, 0).transpose();
    yPt.block<n, 1>(6 + 3 * n, 0) =
        F.block<1, n>(3, 0).transpose(); // dφ2(t,t0) / dt
    yPt.block<n, 1>(6 + 4 * n, 0) = F.block<1, n>(4, 0).transpose();
    yPt.block<n, 1>(6 + 5 * n, 0) = F.block<1, n>(5, 0).transpose();

    return;
}
