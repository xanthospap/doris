#include "astrodynamics.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iersc.hpp"
#include "orbit_integration.hpp"
#include <datetime/dtcalendar.hpp>

constexpr const int Np = 1;

void dso::VariationalEquations(
    double tsec, // seconds from reference epoch (TAI)
    // state and state transition matrix (inertial RF)
    const Eigen::VectorXd &yPhi,
    // state derivative and state transition matrix derivative (inertial RF)
    Eigen::Ref<Eigen::VectorXd> yPhiP,
    // auxiliary parametrs
    dso::IntegrationParameters &params) noexcept {

  // current mjd, TAI
  // const double cmjd = params.mjd_tai + tsec / dso::sec_per_day;
  const dso::TwoPartDate cmjd(params.mjd_tai +
                              dso::TwoPartDate(0e0, tsec / 86400e0));

  // terretrial to celestial for epoch
  Eigen::Matrix<double, 3, 3> rc2i, rpom;
  double era, xlod;
  assert(!gcrs2itrs(cmjd, params.eopLUT, rc2i, era, rpom, xlod));

  // split position and velocity vectors (inertial)
  Eigen::Matrix<double, 3, 1> r = yPhi.block<3, 1>(0, 0);
  Eigen::Matrix<double, 3, 1> v = yPhi.block<3, 1>(3, 0);

  // compute gravity-induced acceleration (we need the position vector in ITRF)
  Eigen::Matrix<double, 3, 3> gpartials;
  Eigen::Matrix<double, 3, 1> r_geo = rcel2ter(r, rc2i, era, rpom);
  Eigen::Matrix<double, 3, 1> gacc;
  test::gravacc3(params.harmonics, r_geo, params.degree, params.harmonics.Re(),
                 params.harmonics.GM(), gacc, gpartials, params.V, params.W);
  // Eigen::Matrix<double, 3, 1> gacc = dso::grav_potential_accel(
  //     r_geo, params.degree, params.order, *(params.Lagrange_V),
  //     *(params.Lagrange_W), params.harmonics, gpartials);

  // fucking crap! gravity acceleration in earth-fixed frame; need to
  // have inertial acceleration!
  gacc = rter2cel(gacc, rc2i, era, rpom);
  const auto rc2ti = Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) * rc2i;
  gpartials = (rpom * rc2ti) * gpartials * (rpom * rc2ti).transpose();

  // third body perturbations, Sun and Moon [m/sec^2] in celestial RF
  Eigen::Matrix<double, 3, 1> rsun; // position of sun, [m] in celestial RF
  Eigen::Matrix<double, 3, 1> rmon; // position of moon, [m] in celestial RF
  Eigen::Matrix<double, 3, 1> sun_acc;
  Eigen::Matrix<double, 3, 1> mon_acc;
  Eigen::Matrix<double, 3, 3> tb_partials;
  dso::SunMoon(cmjd, r, params.GMSun, params.GMMon, sun_acc, mon_acc,
               rsun, rmon, tb_partials);

  // oean tides on geopotential, gravity
  Eigen::Matrix<double, 3, 1> gacc_oc,tmp;
  params.octide->acceleration(cmjd.tai2tt(), /*ut1,*/ r_geo, tmp);
  gacc_oc = rter2cel(tmp, rc2i, era, rpom);
  
  // earth tides on geopotential, gravity
  Eigen::Matrix<double, 3, 1> gacc_ec;
  {
    // Sun and Moon position in ECEF
    const Eigen::Matrix<double, 3, 1> rm_ecef =
        dso::rcel2ter(rmon, rc2i, era, rpom);
    const Eigen::Matrix<double, 3, 1> rs_ecef =
        dso::rcel2ter(rsun, rc2i, era, rpom);
    const auto utc = cmjd.tai2utc();
    dso::EopRecord eops;
    assert(!params.eopLUT.interpolate(utc, eops));
    const auto ut1 =
        dso::TwoPartDate(utc._big, utc._small + eops.dut / 86400e0);
    params.setide->acceleration(cmjd.tai2tt(), ut1, r_geo, rm_ecef, rs_ecef,
                                tmp);
    gacc_ec = rter2cel(tmp, rc2i, era, rpom);
  }

  // Drag
  // Warning only valid for Jason-3
  // get the quaternion
  Eigen::Matrix<double, 3, 1> drag = Eigen::Matrix<double, 3, 1>::Zero();
  Eigen::Matrix<double, 3, 3> ddragdr;
  Eigen::Matrix<double, 3, 3> ddragdv;
  Eigen::Matrix<double, 3, 1> ddragdC;
  Eigen::Quaternion<double> q;
  if (params.qhunt->get_at(cmjd, q)) {
    fprintf(stderr, "ERROR Failed to find quaternion for datetime\n");
    assert(false);
  } else {
    // compute cross-section area
    constexpr const double omegav[] = {0e0, 0e0, iers2010::OmegaEarth};
    const Eigen::Matrix<double, 3, 1> omega{omegav};
    // Velocity relative to the Earth's atmosphere
    const Eigen::Matrix<double, 3, 1> vrel = v - omega.cross(r);
    // normalize
    [[maybe_unused]] const Eigen::Matrix<double, 3, 1> vr = vrel.normalized();
    // loop over flat plates of satellite
    double ProjArea = 0e0;
    for (int i = 0; i < params.numMacroModelComponents; i++) {
      const Eigen::Matrix<double, 3, 1> nb(params.macromodel[i].m_normal);
      const Eigen::Matrix<double, 3, 1> rv = q.conjugate().normalized() * nb;
      const double ctheta = rv.dot(vr);
      if (ctheta > 0e0) {
        ProjArea += params.macromodel[i].m_surf * ctheta;
      }
    }
    // get atmospheric density, using the UTC date
    // TODO for now use TAI date
    assert(
        !params.AtmDataFeed->update_params(cmjd._big, cmjd._small * 86400e0));
    params.AtmDataFeed->set_spatial_from_cartesian(yPhi.block<3, 1>(0, 0));
    dso::nrlmsise00::OutParams aout;
    assert(!params.nrlmsise00->gtd7d(&(params.AtmDataFeed->params_), &aout));
    const double atmdens = aout.d[5];
    Eigen::Matrix<double, 3, 1> drhodr;
    { // approximate arithmetic derivative w.r.t satellite ECEF position
      Eigen::Matrix<double, 3, 1> unitv = Eigen::Matrix<double, 3, 1>::Zero();
      double p1, m1;
      // w.r.t X component
      unitv << 1e0, 0e0, 0e0;
      params.AtmDataFeed->set_spatial_from_cartesian(yPhi.block<3, 1>(0, 0) +
                                                     unitv);
      assert(!params.nrlmsise00->gtd7d(&(params.AtmDataFeed->params_), &aout));
      p1 = aout.d[5];
      unitv << -1e0, 0e0, 0e0;
      params.AtmDataFeed->set_spatial_from_cartesian(yPhi.block<3, 1>(0, 0) +
                                                     unitv);
      assert(!params.nrlmsise00->gtd7d(&(params.AtmDataFeed->params_), &aout));
      m1 = aout.d[5];
      drhodr(0) = ((p1 - atmdens) + (atmdens - m1)) / 2e0;
      // w.r.t Y component
      unitv << 0e0, 1e0, 0e0;
      params.AtmDataFeed->set_spatial_from_cartesian(yPhi.block<3, 1>(0, 0) +
                                                     unitv);
      assert(!params.nrlmsise00->gtd7d(&(params.AtmDataFeed->params_), &aout));
      p1 = aout.d[5];
      unitv << 0e0, -1e0, 0e0;
      params.AtmDataFeed->set_spatial_from_cartesian(yPhi.block<3, 1>(0, 0) +
                                                     unitv);
      assert(!params.nrlmsise00->gtd7d(&(params.AtmDataFeed->params_), &aout));
      m1 = aout.d[5];
      drhodr(1) = ((p1 - atmdens) + (atmdens - m1)) / 2e0;
      // w.r.t Z component
      unitv << 0e0, 0e0, 1e0;
      params.AtmDataFeed->set_spatial_from_cartesian(yPhi.block<3, 1>(0, 0) +
                                                     unitv);
      assert(!params.nrlmsise00->gtd7d(&(params.AtmDataFeed->params_), &aout));
      p1 = aout.d[5];
      unitv << 0e0, 0e0, 1e0;
      params.AtmDataFeed->set_spatial_from_cartesian(yPhi.block<3, 1>(0, 0) +
                                                     unitv);
      assert(!params.nrlmsise00->gtd7d(&(params.AtmDataFeed->params_), &aout));
      m1 = aout.d[5];
      drhodr(2) = ((p1 - atmdens) + (atmdens - m1)) / 2e0;
    }
    drag =
        dso::drag_accel(r, v, ProjArea, *(params.drag_coef), *(params.SatMass),
                        atmdens, drhodr, ddragdr, ddragdv, ddragdC);
  }

  // split state transition and S matrix,
  // yPhi =  | y, F |, size y: 6x1
  //                        F: 6x6
  // but inside the yPhi matrix, they are aranged in a single row, in a
  // column-wise fashion, aka:
  // yPhi = [y_0, y_1, ..., y_5,
  //         F00, F10, ..., F50,
  //         F01, F11, ..., F51,
  //         ...
  //         F05, F15, ..., F55

  // State transition (skip first column which is the state vector)
  // Eigen::Matrix<double, 6, 6> Phi(yPhi.data() + 6);
  Eigen::Matrix<double, 6, 6 + Np> Phi(yPhi.data() + 6);

  // derivative of state transition matrix, aka
  // |   0 (3x3)     I (3x3)   |
  // | da/dr (3x3) da/dv (3x3) |
  // because dv/dr = 0 and
  //         da/dr = I
  Eigen::Matrix<double, 6, 6> dfdy;
  dfdy.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Zero();
  dfdy.block<3, 3>(0, 3) = Eigen::Matrix<double, 3, 3>::Identity();
  dfdy.block<3, 3>(3, 0) = gpartials + tb_partials + ddragdr;
  dfdy.block<3, 3>(3, 3) = /*Eigen::Matrix<double, 3, 3>::Zero()*/ ddragdv;

  // derivative of sensitivity matrix:
  // | 0 (3x3) 0 (3x3) |
  // | 0 (3x3) da/dp   |
  Eigen::Matrix<double, 6, 6 + Np> dfdS =
      Eigen::Matrix<double, 6, 6 + Np>::Zero();
  dfdS.block<3, 1>(3, 6) = ddragdC;

  // Derivative of combined state vector and state transition matrix
  Eigen::Matrix<double, 6, 1 + 6 + Np> yPhip;
  yPhip.block<6, 6 + Np>(0, 1) = dfdy * Phi + dfdS;

  // state derivative (aka [v,a]), in one (first) column
  yPhip.block<3, 1>(0, 0) = v;
  yPhip.block<3, 1>(3, 0) = gacc + sun_acc + mon_acc + drag + gacc_oc + gacc_ec;

  // matrix to vector (column-wise)
  yPhiP = Eigen::VectorXd(
      Eigen::Map<Eigen::VectorXd>(yPhip.data(), yPhip.cols() * yPhip.rows()));

  return;
}
