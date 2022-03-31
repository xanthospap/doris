#ifndef __IMPROVED_RUNGE_KUTTA_NYSTROM_ITERGRATOR_HPP__
#define __IMPROVED_RUNGE_KUTTA_NYSTROM_ITERGRATOR_HPP__

#include "orbit_integrators.hpp"
#include <type_traits>
#include <cstdio>
#ifdef DEBUG
#include <cassert>
#endif

namespace dso {

#if __cplusplus >= 202002L
template <int Order>
class IRKNIntegrator requires orbit_integrators::
    ValidIRKNOrder {
#else
template <int Order,
          typename = std::enable_if_t<
              orbit_integrators::ValidIRKNOrder<Order>::value>>
class IRKNIntegrator {
#endif
private:
  double h;                        ///< step size
  orbit_integrators::ode2_pv *ode; ///< 2nd degree diff. equation

public:
  /// @brief Constructor
  /// @param[in] step_size The step size (aka h)
  /// @param[in] function The 2nd order ODE function; see the 
  ///            orbit_integrators::ode2_pv function alias for more info
  /// @note This function ony assigns; to actually initialize the integrator, 
  /// you should call the initialize() given initial values.
  IRKNIntegrator(double step_size, 
                 orbit_integrators::ode2_pv *function) noexcept
      : h(step_size), ode(function){};

  /// @brief Compute the position and velocity vectors at the next epoch,
  /// aka, t+h
  /// @param[in] t Current epoch
  /// @param[in] r Solution of position vector at t+h
  /// @param[in] r Solution of velocity vector at t+h 
  /// @return t+h
  double step(double t, Vector3 &r, Vector3 &v) noexcept {
  }
}// dso

#endif
