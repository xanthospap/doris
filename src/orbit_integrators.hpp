#ifndef __GENERAL_ORBITIT_ITERGRATOR_HPP__
#define __GENERAL_ORBIT_ITERGRATOR_HPP__

#include "iers2010/matvec.hpp"
#include <type_traits>

namespace dso {

namespace orbit_integrators {
/// @brief Alias for a 2nd order ODE, where the position and velocity vectors
///        are input individually (not as a single state vector). Hence the
///        signature reads:
///        some_function(t, position vec, velocity vec, acceleration vec)
///        and the ODE is used to compute the acceleration given the position 
///        and velocity vectors.
using ode2_pv = void(double, const Vector3 &, const Vector3 &,
                     Vector3 &) noexcept;

/// @brief Valid order for Gauss-Jackson integrators; check at compile-time
#if __cplusplus >= 202002L
    template<int N> concept ValidGausJacksonOrder = N > 1 && N < 12;
#else
    template<int N> struct ValidGaussJacksonOrder {
        static bool const value = N > 1 && N < 12;
    };
#endif

/// @brief Valid order for Improved Runge-Kutta-Nystrom integrators; check at 
/// compile-time
#if __cplusplus >= 202002L
    template<int N> concept ValidIRKNOrder = N > 2 && N < 5;
#else
    template<int N> struct ValidIRKNOrder {
        static bool const value = N > 2 && N < 5;
    };
#endif
} // orbit_integrators
} // dso

#endif
