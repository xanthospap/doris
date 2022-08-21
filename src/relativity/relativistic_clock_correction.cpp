#include "doris_observation_equations.hpp"
#include "geodesy/geodesy.hpp"

/// @brief Compute the relativistic clock correction for a given point 
///        (normally satellite) at r with velocity v (in ECEF).
/// @param[in] recef Position vector of body in ECEF reference frame [m]
/// @param[in] vecef Velocity vector of body in ECEF reference frame [m/sec]
/// @param[in] Re Equatorial radius of the Earth [m]
/// @param[in] GM Gravitational constant [m^3 / s ^2)]
/// @return The quantity U - V^2 / 2 (where U is the potential restricted to
///         the central term, aka μ/r).
/// @see Lemoine etal, 2016, "Precise orbit determination and station
///         position estimation using DORIS RINEX data", section 2.2
///         and 2.5.2
double
dso::relativistic_clock_correction(const Eigen::Matrix<double, 3, 1> &recef,
                                   const Eigen::Matrix<double, 3, 1> &vecef,
                                   double GM) noexcept {
  // potential
  const double U = GM / recef.norm();

  // velocity squared
  const double V2 = vecef.squaredNorm();

#ifdef DEBUG
printf("\t\tDU (beacon) = %.6f\n", U);
#endif

  // return total potential, aka U + V^2 / 2
  return U + V2/2e0;
}

/// @brief Compute the relativistic clock correction for a given point 
///        (normally satellite) at r with velocity v (in ECEF).
/// @param[in] recef Position vector of body in ECEF reference frame [m]
/// @param[in] vecef Velocity vector of body in ECEF reference frame [m/sec]
/// @param[in] Re Equatorial radius of the Earth [m]
/// @param[in] GM Gravitational constant [m^3 / s ^2)]
/// @return The quantity U - V^2 / 2 (where U is the potential including the
///         J2 term).
/// @param[in] J2 Zonal harmonic term J2 in the zero-tide system
/// @see Lemoine etal, 2016, "Precise orbit determination and station
///         position estimation using DORIS RINEX data", section 2.2
///         and 2.5.2
double
dso::relativistic_clock_correction(const Eigen::Matrix<double, 3, 1> &recef,
                                   const Eigen::Matrix<double, 3, 1> &vecef,
                                   double GM, double J2, double Re) noexcept {
  // compute ellipsoidal height
  const auto lfh = dso::car2ell<dso::ellipsoid::grs80>(recef);

  const double latitude = lfh[1]; // [rad]

  // norm of position vector
  const double R = recef.norm();

  // potential, including J2 zonal term (Lemoine et al, Eq. 15)
  const double U =
      (GM / R) *
      (1e0 - (Re / R) * (Re / R) * J2 *
                 (3e0 * std::sin(latitude) * std::sin(latitude) - 1e0) / 2e0);

  // velocity squared
  const double V2 = vecef.squaredNorm();

#ifdef DEBUG
printf("\t\tDU (satellite) = %.6f + %.6f (J2=%.6f)\n", U, V2/2e0, J2);
#endif

  // return total potential, aka U + V^2 / 2
  return U + V2 / 2e0;
}
