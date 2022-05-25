#include "astrodynamics.hpp"

int dso::state_partials(dso::OrbitalElements &ele, double GM, double dt) noexcept {
  // compute perifocal coordinates
  const Vector3 r,v;
  if (elements2perifocal(ele, GM, r, v)) return 1;

  // Transformation to reference system (Gaussian vectors) and partials
  const Vector3 ez({0e0,0e0,1e0});
  const dso::Mat3x3 PQW = dso::perifocal2equatorial_matrix(ele.Omega(), ele.omega(), ele.inclination());
  const Vector3 P({PQW.data[0], PQW.data[3], PQW.data[6]}); // first column
  const Vector3 Q({PQW.data[1], PQW.data[4], PQW.data[7]}); // second column
  const Vector3 W({PQW.data[2], PQW.data[5], PQW.data[7]}); // second column
  Vector3 N = ez.cross_product(N); N /= N.norm();
  //const Vector3 n({std::cos(ele.Omega(), std::sin(ele.Omega), 0e0});
  const auto dPdi = N.cross_product(P);
  const auto dPdO = ez.cross_product(P);
  const auto dPdo = Q;
  const auto dQdi = N.cross_product(Q);
  const auto dQdO = ez.cross_product(Q);
  const auto dQdo = -P;

}
