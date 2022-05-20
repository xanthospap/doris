#include "astrodynamics.hpp"
#include <iers2010/matvec.hpp>

using dso::Vector3;

int
dso::propagate_osculating_elements(double GM, const dso::OrbitalElements &ele_t0,
                                   double dt, dso::OrbitalElements &ele_t) noexcept {
  const double a = ele_t0.semimajor();
  const double e = ele_t0.eccentricity();
  const double i = ele_t0.inclination();
  const double Omega = ele_t0.Omega();
  const double omega = ele_t0.omega();
  const double M0 = ele_t0.mean_anomaly();

  // mean anomaly
  double M;
  if (dt == 0e0) {
    M = M0;
  } else {
    const double n = std::sqrt(GM/(a*a*a));
    M = M0 + n*dt;
  }

  // eccentric anomaly
  int ok;
  const double E = dso::kepler(e, M0, &ok);
  const double sE = std::sin(E);
  const double cE = std::cos(E);

}
