#ifndef __DSO_DORIS_SATELLITE_HPP__
#define __DSO_DORIS_SATELLITE_HPP__

#include "iers2010/matvec.hpp"
#include "iers2010/iersc.hpp"

namespace dso {
  /// @brief Calculate the right ascension α and declination δ from position
  ///        vector.
  /// @see Orbital Mechanics for Engineering Students, H. D. Curtis, (Fourth 
  ///      Edition), Chapter 4.1
  /// @param[in] pos Position vector (x, y, z)
  /// @param[out] a Right_ascension the value of right ascension in the range
  ///               0 to 2π (in [rad])
  /// @param[out] d Declination in the range -π/2 to +π/2 (in [rad])
  int pos2ad(const Vector3 &pos, double &a, double&d) noexcept;

  /// @brief Obtain orbital elements from the state vector.
  /// Applying this algorithm to orbits around other planets or the sun amounts
  /// to defining the frame of reference and substituting the appropriate
  /// gravitational parameter μ(=GM).
  /// This procedure for calculating the orbital elements is not unique.
  /// param[in] pos State vector of in-orbit element (satellite) in [m] and
  ///               [m/s], in the order: (X, Y, Z, Vx, Vy, Vz)
  /// param[out] kepler_ele Kerplerian (orbital) elements, in the following
  ///               order (input array must be at least of size 6):
  ///               [0]: magnitude of (specific) angular momentum [m^2 / s]
  ///               [1]: inclination i, [rad], range: [0, π]
  ///               [2]: right ascension of ascending node Ω, [rad],
  ///                    range [0,2π]
  ///               [3]: magnitude of eccentricity vector, e
  ///               [4]: argument of perigee, ω, [rad], range [0,2π]
  ///               [5]: true anomaly θ, [rad], range [0,2π]
  /// @param[in] GM gravitational parameter [m^3 / s^2]
  /// @note If you want the computation to be performed using [km] instead of
  ///       [m], transform the input quantities accordingly (i.e. state and GM)
  /// @see Orbital Mechanics for Engineering Students, H. D. Curtis, (Fourth
  ///      Edition), Chapter 4
  int state2kepler(const double *state, double *kepler_ele,
                   double GM = iers2010::GMe) noexcept;
  
  /// @brief Given orbital elements compute the state vector in the geocentric 
  //         equatorial frame of reference.
  /// Applying this algorithm to orbits around other planets or the sun amounts 
  /// to defining the frame of reference and substituting the appropriate 
  /// gravitational parameter μ(=GM).
  /// param[in] kepler_ele Kerplerian (orbital) elements, in the following 
  ///               order (input array must be at least of size 6):
  ///               [0]: magnitude of (specific) angular momentum [m^2 / s]
  ///               [1]: inclination i, [rad], range: [0, π]
  ///               [2]: right ascension of ascending node Ω, [rad], 
  ///                    range [0,2π]
  ///               [3]: magnitude of eccentricity vector, e
  ///               [4]: argument of perigee, ω, [rad], range [0,2π]
  ///               [5]: true anomaly θ, [rad], range [0,2π]
  /// param[out] pos State vector of in-orbit elements (satellite) in [m] and 
  ///               [m/s], in the order: (X, Y, Z, Vx, Vy, Vz)
  /// @param[in] GM gravitational parameter [m^3 / s^2]
  /// @note If you want the computation to be performed using [km] instead of 
  ///       [m], transform the input quantities accordingly (i.e. state and GM)
  /// @see Orbital Mechanics for Engineering Students, H. D. Curtis, (Fourth
  ///      Edition), Chapter 4
  int kepler2state(const double *kepler_ele, double *state,
                      double GM = iers2010::GMe) noexcept;

  /// @brief Compute the perturbing acceleration due to the atmospheric drag.
  /// @param[in] r  Satellite position vector in the inertial system [m]
  /// @param[in] v  Satellite velocity vector in the inertial system [m/s]
  /// @param[in] Area Cross-section [m^2]
  /// @param[in] mass Spacecraft mass [kg]
  /// @param[in] CD   Drag coefficient
  /// @param[in] atmdens Atmospheric density
  /// @return Acceleration (a=d^2r/dt^2) [m/s^2]
  /// @see Orbital Mechanics for Engineering Students, H. D. Curtis, (Fourth
  ///      Edition), Chapter 10
  /// @see Satellite Orbits, Models, Methods, Applications, O. Montenbruck, E. 
  ///      Gill, Chapter 3.5
  /// @note Vallado[1] gives a more detailed formula for computing the 
  /// satellite's velocity vector relative to the rotating atmosphere, 
  /// including wind variations. However, such data are usually not available.
  /// [1] Fundamentals of Astrodynamics and Applications, A. Vallado, 
  ///     Chapter 8.6.2
  Vector3 drag_accel(const Vector3 &r, const Vector3 &v, double Area,
                        double mass, double CD, double atmdens) noexcept;
}// dso

#endif
