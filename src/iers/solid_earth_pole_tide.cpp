#include "egravity.hpp"
#include "geodesy/units.hpp"
#include "tides.hpp"
#include <stdexcept>
#ifdef DEBUG
#include <cassert>
#endif

int dso::SolidEarthPoleTide::acceleration(
    const dso::TwoPartDate &mjdtt, double xp_sec, double yp_sec,
    const Eigen::Matrix<double, 3, 1> &rsat, Eigen::Matrix<double, 3, 1> &acc,
    Eigen::Matrix<double, 3, 3> &acc_gradient, double GM, double Re) noexcept 
{
  /* compute disturbing potential, actually only the terms: ΔC21 and ΔS21 */
  const auto dcs = delta_stokes_21(mjdtt, xp_sec, yp_sec);

  /* assign to Stokes coefficients */
  dso::StokesCoeffs dCS(2,2,GM,Re);
  dCS.Cmat().fill_with(0e0);
  dCS.Smat().fill_with(0e0);
  dCS.C(2,1) = dcs.dc21;
  dCS.S(2,1) = dcs.ds21;
  
  /* compute acceleration */
  return dso::gravity_acceleration(dCS, rsat, 2, dCS.Re(), dCS.GM(), acc, acc_gradient);
}
