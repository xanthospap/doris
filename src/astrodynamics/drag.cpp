#include "astrodynamics.hpp"
#include <iers2010/iersc.hpp>

/// @brief Compute acceleration due to atmospheric drag
///
/// @param[in] rsat Satellite position in GCRS/ICRF [m]
/// @param[in] vsat Satellite velocity in GCRS/ICRF [m/sec]
/// @param[in] rbpn GCRs to true-of-date transformation matrix.  
///            The matrix operates in the sense V(date) = rbpn * V(GCRS), 
///            where the vector V(date) is with respect to the true equatorial 
///            triad of date and the vector V(GCRS) is with respect to the 
///            Geocentric Celestial Reference System (IAU, 2000).
///            See iers2010::iau::pnm06
/// @param[in] Area Cross-section [m^2]
/// @param[in] CD   Drag coefficient
/// @param[in] Mass Spacecraft mass [kg]
/// @param[in] Atmospheric density [kg/m^3]
Eigen::Matrix<double, 3, 1>
dso::drag_accel(const Eigen::Matrix<double, 3, 1> &rsat,
                const Eigen::Matrix<double, 3, 1> &vsat,
                const Eigen::Matrix<double, 3, 1> &rbpn, double Area, double CD,
                double Mass, double atmdens) noexcept {
  // earth angular velocity vector [rad/sec]
  constexpr const double omegav[] = {0e0, 0e0, iers2010::OmegaEarth};
  const Eigen::Matrix<double, 3, 1> omega{omegav};

  // position and velocity of satellite to true-of-date system
  const auto rtod = rbpn * rsat;
  const auto vtod = rbpn * vsat;

  // Velocity relative to the Earth's atmosphere
  const auto vrel = vtod - omega.cross(rtod);
  const auto vabs = vrel.norm();

  // Acceleration
  Eigen::Matrix<double, 3, 1> acc =
      -0.5e0 * CD * (Area / Mass) * atmdens * vabs * vrel;
  acc = rbpn.transpose() * acc;

  return acc;
}
