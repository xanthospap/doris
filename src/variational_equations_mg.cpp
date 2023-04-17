#include "astrodynamics.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "orbit_integration.hpp"
#include "iers2010/cel2ter.hpp"
#include <datetime/dtcalendar.hpp>
#include <geodesy/ellipsoid.hpp>

constexpr const int accountforpoletide = true;
constexpr const int m = 2;
constexpr const int n = 6 + m;
/*
 * tsec : seconds from reference time
 *  aka current time t = t0 + tsec
 * 
 * yP0 : Initial condition for state and state transition matrix, in column-
 *  wise oder, i.e.
 *  yP0 = [ y0, Phi(t0,t0) ] 
 *      = [ r0_(3x1), v0_(3x1), 
 *          dr/dr0_(3x3), dr/dv0_(3x3), dr/dp_(3xm)
 *          dv/dr0_(3x3), dv/dv0_(3x3), dv/dp_(3xm)
 *        ]
 *  Initially, the state transition matrix should be set to the identity
 *  matrix, i.e. Phi(t0,t0) = I_(6x6)
 */
namespace {
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
    /*printf("[SHDW ] Changing shadow status from %.3f to %.3f at %.9f\n",
            oc_factor, ShadowFactor, tai.mjd());*/
    oc_factor = ShadowFactor;
  }

  // sun solar radiation pressure, Montenbruck, Sec 3.4
  constexpr const double Ps = 4.56e-6; // Nm^(-2)

  // unit vector directed from spacecrafe to Sun (GCRF)
  const Eigen::Matrix<double, 3, 1> svsun = (sun_gcrf - rgcrf).normalized();

  // get attitude rotation matrix / quaternion
  Eigen::Quaternion<double> q;
  assert(!params.svFrame->get_attitude_quaternion(tai, q));

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

  // scale
  const double rss = (sun_gcrf - rgcrf).norm();
  const double fac = -Cr * ShadowFactor * Ps *
                     std::pow(iers2010::AU / rss, 2e0) / params.svFrame->mass();

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

int drag(const dso::TwoPartDate &tai, const Eigen::Matrix<double, 3, 1> &rgcrf,
         const Eigen::Matrix<double, 3, 1> &vgcrf,
         const Eigen::Matrix<double, 3, 1> &ritrf,
         const dso::Itrs2Gcrs &Rot,
         dso::IntegrationParameters &params, Eigen::Matrix<double, 3, 1> &drag,
         Eigen::Matrix<double, 3, 3> &ddragdr,
         Eigen::Matrix<double, 3, 3> &ddragdv,
         Eigen::Matrix<double, 3, 1> &ddragdC) {
  // get attitude rotation matrix / quaternion
  Eigen::Quaternion<double> q;
  assert(!params.svFrame->get_attitude_quaternion(tai, q));

  // relative velocity
  Eigen::Matrix<double, 3, 1> vrel;
  {
    // Velocity relative to the Earth's atmosphere
    const Eigen::Matrix<double, 3, 1> omega{0e0, 0e0, Rot.omega_earth()};
    vrel = vgcrf - omega.cross(rgcrf);
  }

  // loop over box-wings
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
  S /= params.svFrame->mass();

    // get atmospheric density, using the UTC date
    const dso::TwoPartDate utc = tai.tai2utc();
    const auto lfh = dso::car2ell<dso::ellipsoid::grs80>(ritrf);
    if (params.Dtm20.set(utc, lfh[0], lfh[1], lfh[2]/1e3)) {
      return 1;
    }
    if (params.Dtm20.dtm3()) return 1;
    const double rho = params.Dtm20.total_density_kgm3();
    //printf(" Total Density: %.3e", rho);
    Eigen::Matrix<double, 3, 1> drhodr;
    { // approximate arithmetic derivative w.r.t satellite ECEF position
      const double dr = 5e0;
      // w.r.t X component
      auto dlfh = dso::car2ell<dso::ellipsoid::grs80>(ritrf+dr*Eigen::Matrix<double,3,1>({1e0,0e0,0e0}));
      if (params.Dtm20.set(utc, dlfh[0], dlfh[1], dlfh[2]/1e3)) return 1;
      if (params.Dtm20.dtm3()) return 1;
      double rhop = params.Dtm20.total_density_grcm3();
      dlfh = dso::car2ell<dso::ellipsoid::grs80>(ritrf-dr*Eigen::Matrix<double,3,1>({1e0,0e0,0e0}));
      if (params.Dtm20.set(utc, dlfh[0], dlfh[1], dlfh[2]/1e3)) return 1;
      if (params.Dtm20.dtm3()) return 1;
      double rhom = params.Dtm20.total_density_grcm3();
      drhodr(0) = (rhop-rhom) / (2*dr);
      // w.r.t Y component
      dlfh = dso::car2ell<dso::ellipsoid::grs80>(ritrf+dr*Eigen::Matrix<double,3,1>({0e0,1e0,0e0}));
      if (params.Dtm20.set(utc, dlfh[0], dlfh[1], dlfh[2]/1e3)) return 1;
      if (params.Dtm20.dtm3()) return 1;
      rhop = params.Dtm20.total_density_grcm3();
      dlfh = dso::car2ell<dso::ellipsoid::grs80>(ritrf-dr*Eigen::Matrix<double,3,1>({0e0,1e0,0e0}));
      if (params.Dtm20.set(utc, dlfh[0], dlfh[1], dlfh[2]/1e3)) return 1;
      if (params.Dtm20.dtm3()) return 1;
      rhom = params.Dtm20.total_density_grcm3();
      drhodr(1) = (rhop-rhom) / (2*dr);
      // w.r.t Z component
      dlfh = dso::car2ell<dso::ellipsoid::grs80>(ritrf+dr*Eigen::Matrix<double,3,1>({0e0,0e0,1e0}));
      if (params.Dtm20.set(utc, dlfh[0], dlfh[1], dlfh[2]/1e3)) return 1;
      if (params.Dtm20.dtm3()) return 1;
      rhop = params.Dtm20.total_density_grcm3();
      dlfh = dso::car2ell<dso::ellipsoid::grs80>(ritrf-dr*Eigen::Matrix<double,3,1>({0e0,0e0,1e0}));
      if (params.Dtm20.set(utc, dlfh[0], dlfh[1], dlfh[2]/1e3)) return 1;
      if (params.Dtm20.dtm3()) return 1;
      rhom = params.Dtm20.total_density_grcm3();
      drhodr(2) = (rhop-rhom) / (2*dr);
      // scale (gram/cm3 to kg/m^3)
      drhodr *= 1e-7;
    }
    drag =
        dso::drag_accel(vrel, S, params.Cd, rho, drhodr, ddragdr, ddragdv, ddragdC);
    return 0;
}

}//unnamed namespace

void dso::VariationalEquations_mg(
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

  printf("[ACC] %.9e", cmjd.mjd());

  // split input vector to y, [Φ S]
  Eigen::Matrix<double,6,1> y;
  Eigen::Matrix<double,6,n> FS;
  {
    y = yP0.block<6,1>(0,0).transpose();
    for (int i=0; i<6; i++) {
      FS.row(i) = yP0.block<n,1>(6+n*i,0).transpose();
    }
  }

  // terrestrial to celestial for epoch
  dso::Itrs2Gcrs Rot(cmjd.tai2tt(), &params.eopLUT);

  // split position and velocity vectors GCRF
  Eigen::Matrix<double, 3, 1> r = y.block<3, 1>(0, 0);
  Eigen::Matrix<double, 3, 1> v = y.block<3, 1>(3, 0);

  Eigen::Matrix<double, 3, 1> a  = Eigen::Matrix<double, 3, 1>::Zero();
  Eigen::Matrix<double, 3, 3> dadr = Eigen::Matrix<double, 3, 3>::Zero();
  Eigen::Matrix<double, 3, 3> dadv = Eigen::Matrix<double, 3, 3>::Zero();
  Eigen::Matrix<double, 3, m> dadp = Eigen::Matrix<double, 3, m>::Zero();

  // satellite position at t, in ITRF
  Eigen::Matrix<double, 3, 1> r_itrf = Rot.gcrf2itrf(r);
  { // compute gravity-induced acceleration (we need the position vector in
    // ITRF)
    Eigen::Matrix<double, 3, 3> gpartials;
    Eigen::Matrix<double, 3, 1> gacc;
    test::gravacc3(params.harmonics, r_itrf, params.degree,
                   params.harmonics.Re(), params.harmonics.GM(), gacc,
                   gpartials, params.V, params.W);

    // gravity acceleration in earth-fixed frame; need to have inertial 
    // acceleration!
    gacc = Rot.itrf2gcrf(gacc);
    const auto T = Rot.itrf2gcrf();
    gpartials = T * gpartials * T.transpose();

    a += gacc;
    dadr += gpartials;
    printf(" grav:%.9e", gacc.norm());
    test::gravacc3(params.harmonics, r_itrf, 1,
                   params.harmonics.Re(), params.harmonics.GM(), gacc,
                   gpartials, params.V, params.W);
    printf(" gm:%.9e", Rot.itrf2gcrf(gacc).norm());
  }

  // position of sun/moon, [m] in celestial RF
  Eigen::Matrix<double, 3, 1> rsun,rmon; 
  { // third body perturbations, Sun and Moon [m/sec^2] in celestial RF
    Eigen::Matrix<double, 3, 1> sun_acc;
    Eigen::Matrix<double, 3, 1> mon_acc;
    Eigen::Matrix<double, 3, 3> partials;
    dso::SunMoon(cmjd, r, params.GMSun, params.GMMon, sun_acc, mon_acc, rsun,
                 rmon, partials);
    a += (sun_acc + mon_acc);
    printf(" moon:%.9e sun:%.9e", mon_acc.norm(), sun_acc.norm());
    dadr += partials;
  }

  if (params.setide) { // earth tides on geopotential, gravity
    Eigen::Matrix<double, 3, 1> tacc;
    Eigen::Matrix<double, 3, 3> taccgrad;
    // Sun and Moon position in ECEF
    const Eigen::Matrix<double, 3, 1> rm_ecef = Rot.gcrf2itrf(rmon);
    const Eigen::Matrix<double, 3, 1> rs_ecef = Rot.gcrf2itrf(rsun);
    if (accountforpoletide) {
      params.setide->acceleration(cmjd.tai2tt(), Rot.ut1(),
                                  Rot.eop().xp, Rot.eop().yp, r_itrf, rm_ecef,
                                  rs_ecef, tacc, taccgrad);
    } else {
      params.setide->acceleration(cmjd.tai2tt(), Rot.ut1(), r_itrf, rm_ecef,
                                  rs_ecef, tacc, taccgrad);
    }
    a += Rot.itrf2gcrf(tacc);
    printf(" setide:%.9e", Rot.itrf2gcrf(tacc).norm());
    const auto T = Rot.itrf2gcrf();
    dadr += T * taccgrad * T.transpose();
  }

  if (params.octide) 
  { // oean tides on geopotential, gravity
    Eigen::Matrix<double, 3, 1> tacc;
    params.octide->acceleration(cmjd.tai2tt(), r_itrf, tacc);
    a += Rot.itrf2gcrf(tacc);
    printf(" octide:%.9e", Rot.itrf2gcrf(tacc).norm());
  }

  // drag
  {
    Eigen::Matrix<double, 3, 1> ta;
    Eigen::Matrix<double, 3, 3> tdadr;
    Eigen::Matrix<double, 3, 3> tdadv;
    Eigen::Matrix<double, 3, 1> tdadp;
    if (drag(cmjd.tai2tt(), r, v, r_itrf, Rot, params, ta, tdadr, tdadv, tdadp))
      assert(false);
    a += ta;
    dadr += tdadr;
    dadv += tdadv;
    dadp.block<3,1>(0,0) = tdadp;
    printf(" drag:%.9e", ta.norm());
  }

  {
    Eigen::Matrix<double, 3, 1> ta;
    Eigen::Matrix<double, 3, 3> tdadr;
    Eigen::Matrix<double, 3, 3> tdadv;
    Eigen::Matrix<double, 3, 1> tdadp;
    assert(!solar_radiation_pressure(cmjd, r, rsun, params, params.Cr, ta,
                                     tdadr, tdadv, tdadp));
    a += ta;
    dadr += tdadr;
    dadv += tdadv;
    dadp.block<3,1>(0,1) = tdadp;
    printf(" srp:%.9e", ta.norm());
  }

  // relativity
  {
    Eigen::Matrix<double, 3, 1> ta;
    ta = iers2010::relativistic_correction(r,v);
    a += ta;
    printf(" rel:%.9e", ta.norm());
  }

  Eigen::Matrix<double, 6, 6> F1; // dfdy
  {
    F1.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Zero();
    F1.block<3, 3>(0, 3) = Eigen::Matrix<double, 3, 3>::Identity();
    F1.block<3, 3>(3, 0) = dadr;
    F1.block<3, 3>(3, 3) = dadv;
  }
  
  Eigen::Matrix<double, 6, n> F2; // dfdy
  {
    F2.block<6, 6>(0, 0) = Eigen::Matrix<double, 6, 6>::Zero();
    F2.block<3, m>(0, 6) = Eigen::Matrix<double, 3, m>::Zero();
    F2.block<3, m>(3, 6) = dadp;
  }

  Eigen::Matrix<double, 6, n> VEM = F1 * FS + F2;
  // derivatives
  {
    yPt.block<3,1>(0,0) = v;
    yPt.block<3,1>(3,0) = a;
    for (int i=0; i<6; i++) {
      yPt.block<n,1>(6+i*n,0) = VEM.row(i).transpose();
    }
  }

  printf("\n");
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

void dso::VariationalEquations_ta(
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

  Eigen::Matrix<double, 3, 1> a    = Eigen::Matrix<double, 3, 1>::Zero();
  Eigen::Matrix<double, 3, 3> dadr = Eigen::Matrix<double, 3, 3>::Zero();
  Eigen::Matrix<double, 3, 3> dadv = Eigen::Matrix<double, 3, 3>::Zero();
  Eigen::Matrix<double, 3, m> dadp = Eigen::Matrix<double, 3, m>::Zero();

  // satellite position at t=t0 (ITRF)
  Eigen::Matrix<double, 3, 1> r_itrf = Rot.gcrf2itrf(r);

  { // compute gravity-induced acceleration
    Eigen::Matrix<double, 3, 3> gradient;
    Eigen::Matrix<double, 3, 1> acc;
    test::gravacc3(params.harmonics, r_itrf, params.degree,
                   params.harmonics.Re(), params.harmonics.GM(), acc,
                   gradient, params.V, params.W);

    // transform acceleration and gradient to GCRF
    acc = Rot.itrf2gcrf(acc);
    const auto T = Rot.itrf2gcrf();
    gradient = T * gradient* T.transpose();

    a += acc;
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
    a += (sun_acc + mon_acc);
    dadr += gradient;
  }

  if (params.setide) { 
    // earth tides on geopotential
    Eigen::Matrix<double, 3, 1> acc;
    Eigen::Matrix<double, 3, 3> gradient;
    // Sun and Moon position in ECEF
    const Eigen::Matrix<double, 3, 1> rm_ecef = Rot.gcrf2itrf(rmon);
    const Eigen::Matrix<double, 3, 1> rs_ecef = Rot.gcrf2itrf(rsun);
    if (accountforpoletide) {
      params.setide->acceleration(cmjd.tai2tt(), Rot.ut1(),
                                  Rot.eop().xp, Rot.eop().yp, r_itrf, rm_ecef,
                                  rs_ecef, acc, gradient);
    } else {
      params.setide->acceleration(cmjd.tai2tt(), Rot.ut1(), r_itrf, rm_ecef,
                                  rs_ecef, acc, gradient);
    }
    // acceleration and gradient to GCRF
    a += Rot.itrf2gcrf(acc);
    const auto T = Rot.itrf2gcrf();
    dadr += T * gradient * T.transpose();
  }

  if (params.octide) {
    // oean tides on geopotential, gravity
    Eigen::Matrix<double, 3, 1> acc;
    params.octide->acceleration(cmjd.tai2tt(), r_itrf, acc);
    a += Rot.itrf2gcrf(acc);
  }
  
  // drag
  {
    Eigen::Matrix<double, 3, 1> ta;
    Eigen::Matrix<double, 3, 3> tdadr;
    Eigen::Matrix<double, 3, 3> tdadv;
    Eigen::Matrix<double, 3, 1> tdadp;
    if (drag(cmjd.tai2tt(), r, v, r_itrf, Rot, params, ta, tdadr, tdadv, tdadp))
      assert(false);
    //a += (ta=Eigen::Matrix<double, 3, 1>::Zero());
    //dadr += (tdadr=Eigen::Matrix<double, 3, 3>::Zero());
    //dadv += (tdadv=Eigen::Matrix<double, 3, 3>::Zero());
    //dadp.block<3,1>(0,0) = (tdadp=Eigen::Matrix<double, 3, 1>::Zero());
    a += ta;
    dadr += tdadr;
    dadv += tdadv;
    dadp.block<3,1>(0,0) = tdadp;
    //printf(" %.9f[drag]", Rot.itrf2gcrf(ta).norm());
  }

  {
    Eigen::Matrix<double, 3, 1> ta;
    Eigen::Matrix<double, 3, 3> tdadr;
    Eigen::Matrix<double, 3, 3> tdadv;
    Eigen::Matrix<double, 3, 1> tdadp;
    assert(!solar_radiation_pressure(cmjd, r, rsun, params, params.Cr, ta,
                                     tdadr, tdadv, tdadp));
    a += ta;
    dadr += tdadr;
    dadv += tdadv;
    dadp.block<3,1>(0,1) = tdadp;
    //printf(" %.9f[srp]", ta.norm());
  }
  
  // relativity
  {
    Eigen::Matrix<double, 3, 1> ta;
    ta = iers2010::relativistic_correction(r,v);
    a += ta;
    // printf(" rel:%.9e", ta.norm());
  }
  
  Eigen::Matrix<double,n,n> f123;
  {
    f123.block<3,n>(0,0) = phi1(yP0);
    f123.block<3,n>(3,0) = phi2(yP0);
    if (m) f123.block<m,n>(6,0) = phi3();
  }

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
  Eigen::Matrix<double,n,n> F = A * f123;

  // Assign to yPt
  yPt.block<3,1>(0,0) = v;
  yPt.block<3,1>(3,0) = a;
  yPt.block<n,1>(6+0*n,0) = F.block<1,n>(0,0).transpose(); // dφ1(t,t0) / dt
  yPt.block<n,1>(6+1*n,0) = F.block<1,n>(1,0).transpose();
  yPt.block<n,1>(6+2*n,0) = F.block<1,n>(2,0).transpose();
  yPt.block<n,1>(6+3*n,0) = F.block<1,n>(3,0).transpose(); // dφ2(t,t0) / dt
  yPt.block<n,1>(6+4*n,0) = F.block<1,n>(4,0).transpose();
  yPt.block<n,1>(6+5*n,0) = F.block<1,n>(5,0).transpose();

  return;
}

void dso::noVariationalEquations(
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
    test::gravacc3(params.harmonics, r_geo, params.degree,
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
    if (accountforpoletide) {
      params.setide->acceleration(cmjd.tai2tt(), Rot.ut1(),
                                  Rot.eop().xp, Rot.eop().yp, r_geo, rm_ecef,
                                  rs_ecef, tacc, taccgrad);
    } else {
      params.setide->acceleration(cmjd.tai2tt(), Rot.ut1(), r_geo, rm_ecef,
                                  rs_ecef, tacc, taccgrad);
    }
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
  
  Eigen::Matrix<double, 3, 3> dadv;
  {
    Eigen::Matrix<double, 3, 1> ta;
    Eigen::Matrix<double, 3, 3> tdadr;
    Eigen::Matrix<double, 3, 3> tdadv;
    Eigen::Matrix<double, 3, 1> tdadp;
    assert(!solar_radiation_pressure(cmjd, r, rsun, params, params.Cr, ta,
                                     tdadr, tdadv, tdadp));
    f += ta;
    df += tdadr;
    dadv = tdadv;
    //printf(" %.9f[srp]", ta.norm());
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
