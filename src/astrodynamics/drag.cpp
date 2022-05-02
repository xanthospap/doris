#include "astrodynamics.hpp"

using dso::Vector3;

Vector3 dso::drag_accel(const Vector3 &rsat, const Vector3 &vsat, double Area,
                        double CD, double mass, double atmdens) noexcept {
  // atmospheric "velocity", due to rotation with the earth, Vatm = ω x r
  const Vector3 v_atm{
      {-rsat.y() * iers2010::OmegaEarth, -rsat.x() * iers2010::OmegaEarth}};

  // spacecraft velocity relative to the atmosphere
  const Vector3 v_rel = vsat - v_atm;
  const double v = v_rel.norm();

  // drag force acceleration
  return v_rel * (-0.5e0 * v * atmdens * CD * (Area / mass));
}
