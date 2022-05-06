#ifndef __DSO_ELEMENTARY_ASTRODYNAMICS_HPP__
#define __DSO_ELEMENTARY_ASTRODYNAMICS_HPP__

#include "iers2010/iersc.hpp"
#include "iers2010/matvec.hpp"
#include <limits>

namespace utest {

/// @brief Computes the fractional illumination of a spacecraft in the vicinity
///        of the Earth assuming a cylindrical shadow model
/// @param[in] r_sat Spacecraft position vector [m]
/// @param[in] r_sun Sun position vector [m]
/// @return 0: Spacecraft in Earth shadow; anything other than 0, means that
///            Spacecraft is illuminated by the Sun
/// @ref "Study of satellite shadow function model considering the overlapping
/// parts of Earth shadow and Moon shadow and its application to GPS satellite
/// orbit determination", Zhang et al, 2018
double conical_shadow(const dso::Vector3 &r_sat,
                      const dso::Vector3 &r_sun) noexcept;

/// @brief Computes the fractional illumination of a spacecraft in the vicinity
///        of the Earth assuming a cylindrical shadow model
/// @param[in] r_sat Spacecraft position vector [m]
/// @param[in] r_sun Sun position vector [m]
/// @return 0: Spacecraft in Earth shadow; anything other than 0, means that
///            Spacecraft is illuminated by the Sun
double vallado_shadow(const dso::Vector3 &r_sat,
                      const dso::Vector3 &r_sun) noexcept;

/// @brief Computes the fractional illumination of a spacecraft in the vicinity
///        of the Earth assuming a cylindrical shadow model
/// @param[in] r_sat Spacecraft position vector [m]
/// @param[in] r_sun Sun position vector [m]
/// @return 0: Spacecraft in Earth shadow; anything other than 0, means that
///            Spacecraft is illuminated by the Sun
/// @note Only has two states, 1 or 0
/// @ref Satellite
double montebruck_shadow(const dso::Vector3 &r_sat,
                         const dso::Vector3 &r_sun) noexcept;

double bernese_shadow(const dso::Vector3 &r_sat,
                      const dso::Vector3 &r_sun) noexcept;
double bernese_shadow1(const dso::Vector3 &r_sat,
                       const dso::Vector3 &r_sun) noexcept;
} // namespace utest

namespace dso {

struct OrbitalElements {
  /// [0]: Semimajor axis a
  /// [1]: Eccentricity e
  /// [2]: Inclination i, [rad], range: [0, π]
  /// [3]: Longitude of ascending node Ω, [rad], range [0,2π]
  /// [4]: Argument of pericenter, ω, [rad], range [0,2π]
  /// [5]: Mean Anomaly M [rad]
  /// [6]: True Anomaly θ [rad], range [0,2π]
  double elements[7];
  constexpr double& semimajor() noexcept {return elements[0];}
  constexpr double& eccentricity() noexcept {return elements[1];}
  constexpr double& inclination() noexcept {return elements[2];}
  constexpr double& Omerga() noexcept {return elements[3];}
  constexpr double& omega() noexcept {return elements[4];}
  constexpr double& mean_anomaly() noexcept {return elements[5];}
  constexpr double& true_anomaly() noexcept {return elements[6];}
  constexpr double semimajor()    const noexcept {return elements[0];}
  constexpr double eccentricity() const noexcept {return elements[1];}
  constexpr double inclination()  const noexcept {return elements[2];}
  constexpr double Omerga()       const noexcept {return elements[3];}
  constexpr double omega()        const noexcept {return elements[4];}
  constexpr double mean_anomaly() const noexcept {return elements[5];}
  constexpr double true_anomaly() const noexcept {return elements[6];}
};
int state2elements(const dso::Vector3 &r, const dso::Vector3 &v,
                   OrbitalElements &elements, double GM) noexcept;

/// @brief Solve Kepler's equation iteratively via Newton's method.
/// @param[in] e Orbit eccentricity
/// @param[in] M mean anomaly [radians]
/// @param[in] tolerance_rad End iterations when: f(E) < tolerance_rad
///            units in [radians] (Note that f(E) = E - esinE - M)
/// @param[out] ok Anything other than 0, denotes a convergence error and
///            the result should not be used.
/// @return Eccentric Anomaly, E [radians]
/// @note For small eccentricities, the starting value for the iteration is
///       E0 = M, since E only differs from M by a term of order e. For
///       highly eccentric orbits however (e.g. e > 0.8), the iteration
///       should start with E0 = π to avoid convergence problems.
/// @see Montenbruck et al, 2000, Eq. 2.42
double
kepler(double e, double M, int &ok,
       double tolerance_rad = 1e2 *
                              std::numeric_limits<double>::epsilon()) noexcept;

/// @brief Solve Kepler's equation iteratively via Newton's method.
/// @param[in] e Orbit eccentricity
/// @param[in] M mean anomaly [radians]
/// @param[in] tolerance_rad End iterations when:
///            |E_n+1 - E_n| < tolerance_rad
///            units in [radians
/// @param[out] ok Anything other than 0, denotes a convergence error and
///            the result should not be used.
/// @return Eccentric Anomaly, E [radians]
/// @note Starting value is:
///      -π < M < 0 or M > π => E = M - e
///           else           => E = M + e
/// @warning Only works for Elliptical orbits; do not use for Parabolic and
///          Hyperbolic orbits.
/// @see  Vallado, 2.2.5
double kepler_vallado(
    double e, double M, int &ok,
    double tolerance_rad = 1e2 *
                           std::numeric_limits<double>::epsilon()) noexcept;

/// @brief Calculate the right ascension α and declination δ from position
///        vector. The position vector must be given in the (non-rotating) 
///        geocentric equatorial frame.
/// @see Orbital Mechanics for Engineering Students, H. D. Curtis, (Fourth
///      Edition), Chapter 4.1
/// @param[in] pos Position vector (x, y, z) in the geocentric equatorial 
///               frame
/// @param[out] a Right_ascension the value of right ascension in the range
///               0 to 2π (in [rad])
/// @param[out] d Declination in the range -π/2 to +π/2 (in [rad])
int pos2ad(const Vector3 &pos, double &a, double &d) noexcept;

/// @brief Obtain orbital elements from the state vector.
/// The state vector (position and velocity) must be given in the geocentric
/// equatorial frame.
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
namespace alternatives {
int state2kepler_vallado(const double *state, double *kepler_ele,
                         double GM = iers2010::GMe) noexcept;
int state2kepler_montenbruck(const double *state, double *kepler_ele,
                         double GM = iers2010::GMe) noexcept;
}

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
Vector3 drag_accel(const Vector3 &r, const Vector3 &v, double Area, double mass,
                   double CD, double atmdens) noexcept;
} // namespace dso

#endif
