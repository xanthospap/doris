#include "astrodynamics.hpp"
#include "geodesy/units.hpp"

/* dYdA */
Eigen::Matrix<double, 6, 6> dso::state_partials(double GM,
                                                const dso::OrbitalElements &ele,
                                                double dt) noexcept {
  // compute perifocal coordinates
  Eigen::Matrix<double, 6, 1> state = dso::elements2perifocal(GM,ele,dt);
  const double x = state(0);
  const double y = state(1);
  const double vx = state(3);
  const double vy = state(4);
  const double r = state.block<3, 1>(0, 0).norm();

  // Transformation to reference system (Gaussian vectors) and partials
  Eigen::Matrix3d PQW(
      Eigen::AngleAxisd(ele.Omega(), Eigen::Vector3d::UnitZ()) *
      Eigen::AngleAxisd(ele.inclination(), Eigen::Vector3d::UnitX()) *
      Eigen::AngleAxisd(ele.omega(), Eigen::Vector3d::UnitZ()));
  const auto P = PQW.block<3, 1>(0, 0); // first column
  const auto Q = PQW.block<3, 1>(0, 1); // second column
  const auto W = PQW.block<3, 1>(0, 2); // third column
  
  auto N = Eigen::Vector3d::UnitZ().cross(W);
  N = N / N.norm();

  // const Vector3 n({std::cos(ele.Omega(), std::sin(ele.Omega), 0e0});
  const auto dPdi = N.cross(P);
  const auto dPdO = Eigen::Vector3d::UnitZ().cross(P);
  const auto dPdo = Q;
  const auto dQdi = N.cross(Q);
  const auto dQdO = Eigen::Vector3d::UnitZ().cross(Q);
  const auto dQdo = -P;

  // Partials w.r.t. inlcination, node and argument of pericenter
  Eigen::Matrix<double, 6, 1> dYdi;
  dYdi << x * dPdi + y * dQdi, vx * dPdi + vy * dQdi;
  Eigen::Matrix<double, 6, 1> dYdO;
  dYdO << x * dPdO + y * dQdO, vx * dPdO + vy * dQdO;
  Eigen::Matrix<double, 6, 1> dYdo;
  dYdo << x * dPdo + y * dQdo, vx * dPdo + vy * dQdo;

  // Partials w.r.t. semimajor axis, eccentricity and mean anomaly at time dt
  const double a = ele.semimajor();
  const double e = ele.eccentricity();
  const double fac = std::sqrt((1e0 - e) * (1e0 + e));
  const double n = std::sqrt(GM / (a * a * a));
  
  Eigen::Matrix<double, 6, 1> dYda;
  dYda << ((x / a) * P + (y / a) * Q),
      ((-vx / (2e0 * a)) * P + (-vy / (2e0 * a)) * Q);
  
  Eigen::Matrix<double, 6, 1> dYde;
  dYde << (-a - std::pow(y / fac, 2) / r) * P + (x * y / (r * fac * fac)) * Q,
      (vx * (2 * a * x + e * std::pow(y / fac, 2)) / (r * r)) * P +
          ((n / fac) * std::pow(a / r, 2) *
           (x * x / r - std::pow(y / fac, 2) / a)) *
              Q;
  
  Eigen::Matrix<double, 6, 1> dYdM;
  dYdM << ((vx * P + vy * Q) / n),
      ((-n * std::pow(a / r, 3e0)) * (x * P + y * Q));

  // Derivative of mean anomaly at time dt w.r.t. the semimajor axis at epoch
  const double dMda = -1.5e0 * (n / a) * dt;

  // Combined partial derivative matrix of state with respect to epoch elements
  Eigen::Matrix<double, 6, 6> dYdA;
  for (int k = 0; k < 6; k++) {
    dYdA(k, 0) = dYda(k) + dYdM(k) * dMda;
    dYdA(k, 1) = dYde(k);
    dYdA(k, 2) = dYdi(k);
    dYdA(k, 3) = dYdO(k);
    dYdA(k, 4) = dYdo(k);
    dYdA(k, 5) = dYdM(k);
  }

  return dYdA;
}
