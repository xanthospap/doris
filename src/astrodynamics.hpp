#ifndef __DSO_ELEMENTARY_ASTRODYNAMICS_HPP__
#define __DSO_ELEMENTARY_ASTRODYNAMICS_HPP__

#include <limits>

namespace dso {
  /// @brief Solve Kepler's equation iteratively via Newton's method.
  /// @param[in] e Orbit eccentricity
  /// @param[in] M mean anomaly [radians]
  /// @param[in] tolerance_rad End iterations when: f(E) < tolerance_rad
  ///            units in [radians] (Note that f(E) = E - esinE - M)
  /// @return Eccentric Anomaly, E [radians]
  /// @note For small eccentricities, the starting value for the iteration is
  ///       E0 = M, since E only differs from M by a term of order e. For 
  ///       highly eccentric orbits however (e.g. e > 0.8), the iteration 
  ///       should start with E0 = π to avoid convergence problems.
  /// @see Montenbruck et al, 2000, Eq. 2.42
double
kepler(double e, double M,
       double tolerance_rad = 1e2 *
                              std::numeric_limits<double>::epsilon()) noexcept;
}// dso

#endif
