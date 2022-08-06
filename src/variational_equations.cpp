#include "orbit_integration.hpp"

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
  Eigen::Matrix<double, 3, 3> dt2c;
  Eigen::Matrix<double, 3, 3> t2c(
      dso::itrs2gcrs(cmjd, params.eopLUT, dt2c));

  // split position and velocity vectors (inertial)
  Eigen::Matrix<double, 3, 1> r = yPhi.block<3, 1>(0, 0);
  Eigen::Matrix<double, 3, 1> v = yPhi.block<3, 1>(3, 0);

  // compute gravity-induced acceleration
  Eigen::Matrix<double, 3, 3> gpartials;
  Eigen::Matrix<double, 3, 1> r_geo = t2c.transpose() * r;
  Eigen::Matrix<double, 3, 1> gacc = dso::grav_potential_accel(
      r_geo, params.degree, params.order, *(params.Lagrange_V), *(params.Lagrange_W),
      params.harmonics, gpartials);

  // fucking crap! gravity acceleration in earth-fixed frame; need to
  // have inertial acceleration!
  gacc = t2c * gacc;
  gpartials = t2c.transpose() * gpartials * t2c;

  // third body perturbations, Sun and Moon [m/sec^2] in celestial RF
  Eigen::Matrix<double, 3, 1> rsun; // position of sun, [m] in celestial RF
  Eigen::Matrix<double, 3, 1> sun_acc;
  Eigen::Matrix<double, 3, 1> mon_acc;
  Eigen::Matrix<double, 3, 3> tb_partials;
  SunMoon(cmjd, r, sun_acc, mon_acc, rsun,
          tb_partials);

  // Drag
  // set input parameters (time and spatial)
  //Eigen::Matrix<double, 3, 1> drag_acc = Eigen::Matrix<double, 3, 1>::Zero();
  //if (include_drag) {
  //  int msise_mjd;
  //  double msise_secday;
  //  setNrlmsise00Params(cmjd, r_geo, msise_mjd, msise_secday, params->msise_in);
  //  //printf("\tparams set for setNrlmsise00Params ...\n");
  //  // update flux data
  //  if (params->msise_in->update_params(msise_mjd, msise_secday)) {
  //    fprintf(stderr, "ERROR. Failed to update flux/Ap data from drag!\n");
  //  }
  //  //printf("\tparams updated in hunter ...\n");
  //  // get density (kg/m^3)
  //  double density;
  //  {
  //    dso::nrlmsise00::OutParams out;
  //    //params->msise_in->params_.dump_params();
  //    params->msise->gtd7d(&params->msise_in->params_, &out);
  //    density = out.d[5];
  //    //printf("\tdensity value is: %+.15e\n", out.d[5]);
  //  }
  //  // we also need to true-to-date matrix (for relavite velocity)
  //  Eigen::Matrix<double, 3, 3> r_tof;
  //  {
  //    double mjd_days;
  //    const double taif = std::modf(cmjd, &mjd_days);
  //    const double ttf = taif + (32184e-3 / 86400e0);
  //    const double mjd_tt = mjd_days + ttf;
  //    const auto rnpb = iers2010::sofa::pnm06a(dso::mjd0_jd, mjd_tt);
  //    Eigen::Matrix<double, 3, 3> r_toft(rnpb.data);
  //    r_tof = r_toft.transpose();
  //  }
  //  const auto drag =
  //      dso::drag_accel(r, v, r_tof, Area_H2C, CD, Mass_H2C, density);
  //  drag_acc = drag;
  //  //printf("\tdrag computed, done!\n");
  //}

  // SRP
  //Eigen::Matrix<double, 3, 1> srp = Eigen::Matrix<double, 3, 1>::Zero();
  //if (include_srp) {
  //  dso::Vector3 rV({r(0), r(1), r(2)});
  //  dso::Vector3 sV({rsun(0), rsun(1), rsun(2)});
  //  if (utest::montebruck_shadow(rV, sV)) {
  //    srp = dso::solar_radiation_acceleration(r, rsun, Area_H2C, Mass_H2C,
  //                                            Cr_H2C);
  //  }
  //}

  // State transition (skip first column which is the state vector)
  Eigen::Matrix<double, 6, 6> Phi(yPhi.data() + 6);

  // derivative of state transition matrix, aka
  // | v (3x3)   0 (3x3)     I (3x3)   |
  // | a (3x1) da/dr (3x3) da/dv (3x3) |
  // because dv/dr = 0 and
  //         da/dr = I
  Eigen::Matrix<double, 6, 6> dfdy;
  dfdy.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Zero();
  dfdy.block<3, 3>(0, 3) = Eigen::Matrix<double, 3, 3>::Identity();
  dfdy.block<3, 3>(3, 0) = gpartials + tb_partials;
  dfdy.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Zero();

  // Derivative of combined state vector and state transition matrix
  // dPhi = dfdy * Phi;
  Eigen::Matrix<double, 6, 7> yPhip;
  yPhip.block<6, 6>(0, 1) = dfdy * Phi;

  // state derivative (aka [v,a]), in one (first) column
  yPhip.block<3, 1>(0, 0) = v;
  yPhip.block<3, 1>(3, 0) = gacc + sun_acc + mon_acc;

  // matrix to vector (column-wise)
  yPhiP = Eigen::VectorXd(
      Eigen::Map<Eigen::VectorXd>(yPhip.data(), yPhip.cols() * yPhip.rows()));

  //printf("Accelerations: \n");
  //printf("Gravity     Sun          Moon         srp          drag\n");
  //for (int i=0; i<3; i++) {
  //  printf("%+12.9f %+12.9f %+12.9f %+12.9f %+12.9f\n", gacc(i), sun_acc(i), mon_acc(i), srp(i), drag_acc(i));
  //}

  return;
}
