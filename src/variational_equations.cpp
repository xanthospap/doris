#include "orbit_integration.hpp"
#include "iers2010/iersc.hpp"
#include "astrodynamics.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"

constexpr const int Np = 1;

void dso::VariationalEquations(
    double tsec, // TAI
    // state and state transition matrix (inertial RF)
    const Eigen::VectorXd &yPhi,
    // state derivative and state transition matrix derivative (inertial RF)
    Eigen::Ref<Eigen::VectorXd> yPhiP,
    // auxiliary parametrs
    dso::IntegrationParameters &params) noexcept {

  // current mjd, TAI
  const double cmjd = params.mjd_tai + tsec / dso::sec_per_day;

  // terretrial to celestial for epoch
#ifdef ABCD
  Eigen::Matrix<double, 3, 3> dt2c;
  Eigen::Matrix<double, 3, 3> t2c(dso::itrs2gcrs(cmjd, params.eopLUT, dt2c));
#else
  Eigen::Matrix<double, 3, 3> rc2i, rpom;
  double era;
  assert(!gcrs2itrs(cmjd, params.eopLUT, rc2i, era, rpom));
  // const auto t2c = (rpom * rc2ti).transpose() ;
#endif
  //{
  //  printf("Terrestrial to Celestail Matrix:\n");
  //  for (int i=0; i<3; i++) {
  //    for (int j=0; j<3; j++) {
  //      printf(" %+.6f ", t2c(i,j));
  //    }
  //    printf("\n");
  //  }
  //}

  // split position and velocity vectors (inertial)
  Eigen::Matrix<double, 3, 1> r = yPhi.block<3, 1>(0, 0);
  Eigen::Matrix<double, 3, 1> v = yPhi.block<3, 1>(3, 0);

  // compute gravity-induced acceleration (we need the position vector in ITRF)
  Eigen::Matrix<double, 3, 3> gpartials;
#ifdef ABCD
  Eigen::Matrix<double, 3, 1> r_geo = t2c.transpose() * r;
#else
  Eigen::Matrix<double, 3, 1> r_geo = rcel2ter(r, rc2i, era, rpom);
#endif
  Eigen::Matrix<double, 3, 1> gacc = dso::grav_potential_accel(
      r_geo, params.degree, params.order, *(params.Lagrange_V),
      *(params.Lagrange_W), params.harmonics, gpartials);

  // fucking crap! gravity acceleration in earth-fixed frame; need to
  // have inertial acceleration!
  //printf(">> ITRF acc: %+.9f %+.9f %+.9f\n", gacc(0), gacc(1), gacc(2));
#ifdef ABCD
  gacc = t2c * gacc;
  gpartials = t2c.transpose() * gpartials * t2c;
  //for (int i=0; i<3; i++) {
  //  for (int j=0; j<3; j++) {
  //    printf(" %+.6f ", t2c.transpose()(i,j));
  //  }
  //  printf("\t\t");
  //  for (int j=0; j<3; j++) {
  //    printf(" %+.6f ",t2c(i,j));
  //  }
  //  printf("\n");
  //}
#else
  gacc = rter2cel(gacc, rc2i, era, rpom);
  const auto rc2ti = Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) * rc2i;
  gpartials = (rpom * rc2ti) * gpartials * (rpom * rc2ti).transpose();
  //for (int i=0; i<3; i++) {
  //  for (int j=0; j<3; j++) {
  //    printf(" %+.6f ",(rpom * rc2ti)(i,j));
  //  }
  //  printf("\t\t");
  //  for (int j=0; j<3; j++) {
  //    printf(" %+.6f ",(rpom * rc2ti).transpose()(i,j));
  //  }
  //  printf("\n");
  //}
#endif
  //printf(">> GCRF acc: %+.9f %+.9f %+.9f\n", gacc(0), gacc(1), gacc(2));

  // third body perturbations, Sun and Moon [m/sec^2] in celestial RF
  Eigen::Matrix<double, 3, 1> rsun; // position of sun, [m] in celestial RF
  Eigen::Matrix<double, 3, 1> sun_acc;
  Eigen::Matrix<double, 3, 1> mon_acc;
  Eigen::Matrix<double, 3, 3> tb_partials;
  dso::SunMoon(cmjd, r, params.GMSun, params.GMMon, sun_acc, mon_acc, rsun,
               tb_partials);

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
    [[maybe_unused]]const Eigen::Matrix<double, 3, 1> vr = vrel.normalized();
    // loop over flat plates of satellite
    double ProjArea = 0e0;
    for (int i=0; i<params.numMacroModelComponents; i++) {
      const Eigen::Matrix<double, 3, 1> nb(params.macromodel[i].m_normal);
      const Eigen::Matrix<double, 3, 1> rv = q.conjugate().normalized() * nb;
      const double ctheta = rv.dot(vr);
      if (ctheta > 0e0) {
        ProjArea += params.macromodel[i].m_surf * ctheta;
      }
    }
    // get atmospheric density, using the UTC date
    double imjd; 
    int eid;
    double utc_fday = std::modf(cmjd, &imjd);
    const int leap_sec = dso::dat(cmjd, eid);
    utc_fday -= leap_sec / (86400e0 + (double)eid);
    if (utc_fday > 1e0) {
      utc_fday -= 1e0;
      ++imjd;
    } else if (utc_fday < 1e0) {
      utc_fday += 1e0;
      --imjd;
    }
    assert(!params.AtmDataFeed->update_params(imjd, utc_fday*86400e0));
    params.AtmDataFeed->set_spatial_from_cartesian(yPhi.block<3, 1>(0, 0));
    dso::nrlmsise00::OutParams aout;
    assert(!params.nrlmsise00->gtd7d(&(params.AtmDataFeed->params_), &aout));
    const double atmdens = aout.d[5];
    Eigen::Matrix<double,3,1> drhodr;
    { // approximate arithmetic derivative w.r.t satellite ECEF position
      Eigen::Matrix<double,3,1> unitv = Eigen::Matrix<double,3,1>::Zero();
      double p1,m1;
      // w.r.t X component
      unitv << 1e0, 0e0, 0e0;
      params.AtmDataFeed->set_spatial_from_cartesian(yPhi.block<3, 1>(0, 0) + unitv);
      assert(!params.nrlmsise00->gtd7d(&(params.AtmDataFeed->params_), &aout));
      p1 = aout.d[5];
      unitv << -1e0, 0e0, 0e0;
      params.AtmDataFeed->set_spatial_from_cartesian(yPhi.block<3, 1>(0, 0) + unitv);
      assert(!params.nrlmsise00->gtd7d(&(params.AtmDataFeed->params_), &aout));
      m1 = aout.d[5];
      drhodr(0) = ((p1-atmdens) + (atmdens-m1)) / 2e0;
      // w.r.t Y component
      unitv << 0e0, 1e0, 0e0;
      params.AtmDataFeed->set_spatial_from_cartesian(yPhi.block<3, 1>(0, 0) + unitv);
      assert(!params.nrlmsise00->gtd7d(&(params.AtmDataFeed->params_), &aout));
      p1 = aout.d[5];
      unitv << 0e0, -1e0, 0e0;
      params.AtmDataFeed->set_spatial_from_cartesian(yPhi.block<3, 1>(0, 0) + unitv);
      assert(!params.nrlmsise00->gtd7d(&(params.AtmDataFeed->params_), &aout));
      m1 = aout.d[5];
      drhodr(1) = ((p1-atmdens) + (atmdens-m1)) / 2e0;
      // w.r.t Z component
      unitv << 0e0, 0e0, 1e0;
      params.AtmDataFeed->set_spatial_from_cartesian(yPhi.block<3, 1>(0, 0) + unitv);
      assert(!params.nrlmsise00->gtd7d(&(params.AtmDataFeed->params_), &aout));
      p1 = aout.d[5];
      unitv << 0e0, 0e0, 1e0;
      params.AtmDataFeed->set_spatial_from_cartesian(yPhi.block<3, 1>(0, 0) + unitv);
      assert(!params.nrlmsise00->gtd7d(&(params.AtmDataFeed->params_), &aout));
      m1 = aout.d[5];
      drhodr(2) = ((p1-atmdens) + (atmdens-m1)) / 2e0;
    }
  //  //printf("Note: Using drag coefficient=%.3f\n", params.get_drag_coefficient());
    //drag = dso::drag_accel(r, v, ProjArea, /*params.get_drag_coefficient()*/2e0,
    //                       *(params.SatMass), atmdens/*, drhodr, ddragdr, ddragdv,
    //                       ddragdC*/);
    drag = dso::drag_accel(r, v, ProjArea, *(params.drag_coef),
                           *(params.SatMass), atmdens, drhodr, ddragdr, ddragdv,
                           ddragdC);
  }

  // SRP
  // Eigen::Matrix<double, 3, 1> srp = Eigen::Matrix<double, 3, 1>::Zero();
  // if (include_srp) {
  //  dso::Vector3 rV({r(0), r(1), r(2)});
  //  dso::Vector3 sV({rsun(0), rsun(1), rsun(2)});
  //  if (utest::montebruck_shadow(rV, sV)) {
  //    srp = dso::solar_radiation_acceleration(r, rsun, Area_H2C, Mass_H2C,
  //                                            Cr_H2C);
  //  }
  //}
  //

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
  Eigen::Matrix<double, 6, 6+Np> Phi(yPhi.data() + 6);

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
  Eigen::Matrix<double, 6, 6+Np> dfdS = Eigen::Matrix<double, 6, 6+Np>::Zero();
  dfdS.block<3, 1>(3, 6) = ddragdC;

  // Derivative of combined state vector and state transition matrix
  Eigen::Matrix<double, 6, 1+6+Np> yPhip;
  yPhip.block<6,6+Np>(0,1) = dfdy * Phi + dfdS;
  
  // state derivative (aka [v,a]), in one (first) column
  yPhip.block<3, 1>(0, 0) = v;
  yPhip.block<3, 1>(3, 0) = gacc + sun_acc + mon_acc + drag;

  // matrix to vector (column-wise)
  yPhiP = Eigen::VectorXd(
      Eigen::Map<Eigen::VectorXd>(yPhip.data(), yPhip.cols() * yPhip.rows()));

  return;
}
