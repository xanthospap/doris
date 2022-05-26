#include "astrodynamics.hpp"

int dso::state_partials(dso::OrbitalElements &ele, double GM,
                        double dt) noexcept {
  // compute perifocal coordinates
  const Vector3 r, v;
  if (elements2perifocal(ele, GM, r, v))
    return 1;

  // Transformation to reference system (Gaussian vectors) and partials
  const Vector3 ez({0e0, 0e0, 1e0});
  const dso::Mat3x3 PQW = dso::perifocal2equatorial_matrix(
      ele.Omega(), ele.omega(), ele.inclination());
  const Vector3 P({PQW.data[0], PQW.data[3], PQW.data[6]}); // first column
  const Vector3 Q({PQW.data[1], PQW.data[4], PQW.data[7]}); // second column
  const Vector3 W({PQW.data[2], PQW.data[5], PQW.data[7]}); // second column
  Vector3 N = ez.cross_product(N);
  N /= N.norm();
  // const Vector3 n({std::cos(ele.Omega(), std::sin(ele.Omega), 0e0});
  const auto dPdi = N.cross_product(P);
  const auto dPdO = ez.cross_product(P);
  const auto dPdo = Q;
  const auto dQdi = N.cross_product(Q);
  const auto dQdO = ez.cross_product(Q);
  const auto dQdo = -P;

  // perifocal coordinates
  Vector3 r_perif, v_perif;
  if (dso::elements2perifocal(ele, GM, r_perif, v_perif))
    return 1;
  
  // Partials w.r.t. inlcination, node and argument of pericenter
  const auto drdi = r_perif.x()*dPdi + r_perif.y()*dQdi;
  const auto dvdi = v_perif.x()*dPdi + v_perif.y()*dQdi;
  const auto drdO = r_perif.x()*dPdO + r_perif.y()*dQdO;
  const auto dvdO = v_perif.x()*dPdO + v_perif.y()*dQdO;
  const auto drdo = r_perif.x()*dPdo + r_perif.y()*dQdo;
  const auto dvdo = v_perif.x()*dPdo + v_perif.y()*dQdo;

  // Partials w.r.t. semimajor axis, eccentricity and mean anomaly at time dt
  const double r = r_perif.norm();
  const double a = ele.semimajor();
  const double e = ele.eccentricity();
  const double fac = (1e0 - e) * (1e0 + e);
  const auto drda = (r_perif.x() / a) * P + (r_perif.y() / a) * Q;
  const auto dvda =
      (-v_perif.x() / (2e0 * a)) * P + (-v_perif.y() / (2e0 * a)) * Q;
  const auto drde = (-a - (r_perif.y / fac) * (r_perif.y / fac) / r) * P +
                    (r_perif.x() * r_perif.y() / (r * fac * fac)) * Q;
  const double n = std::sqrt(GM / (a * a * a));
  const auto dvde =
      (v_perif.x() *
       (2e0 * a * r_perif.x() + e * std::pow(r_perif.y() / fac, 2e0)) /
       (r * r)) *
          P +
      ((n / fac) * std::pow(a / r, 2e0) *
       (r_perif.x() * r_perif.x() / r - std::pow(r_perif.y() / fac, 2e0) / a)) *
          Q;
  const auto drdM = (v_perif.x()*P+v_perif.y()*Q)/n;
  const auto dvdM = (-n*std::pow(a/r,3e0))*(r_perif.x()*P + r_perif.y()*Q);

  // Derivative of mean anomaly at time dt w.r.t. the semimajor axis at epoch
  const double dMda = -1.5e0 * (n/a) * dt;

  // Combined partial derivative matrix of state with respect to epoch elements
  Eigen::Matrix<double, 6, 6> dYdA;
  for (int i=0; i<3; i++) {
    dYdA(i,0) = drda.data[j] + drdM.data[j] * dMda;
    dYdA(i,1) = drde.data[j];
    dYdA(i,2) = drdi.data[j];
  for (int i=3; i<6; i++) {
    dYdA(i,3) = dvda.data[j] + dvdM.data[j] * dMda;
    dYdA(i,4) = dvde.data[j];
    dYdA(i,5) = dvdi.data[j];
  }

  return dYdA;
}
