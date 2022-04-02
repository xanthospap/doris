#ifndef __IMPROVED_RUNGE_KUTTA_NYSTROM_ITERGRATOR_HPP__
#define __IMPROVED_RUNGE_KUTTA_NYSTROM_ITERGRATOR_HPP__

#include "orbit_integrators.hpp"
#include "orbit_integrators_coefficients.hpp"
#include <type_traits>
#include <cstdio>
#ifdef DEBUG
#include <cassert>
#endif

namespace dso {

template<int N, int M>
class RungeKuttaNystromIntegrator {
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
  RungeKuttaNystromIntegrator(double step_size, 
                 orbit_integrators::ode2_pv *function) noexcept
      : h(step_size), ode(function){};

  /// @brief Compute the position and velocity vectors at the next epoch,
  /// aka, t+h
  /// @param[in] t Current epoch
  /// @param[in] r Solution of position vector at t+h
  /// @param[in] r Solution of velocity vector at t+h 
  /// @return t+h
  double step(double t, Vector3 &r, Vector3 &v) noexcept {
    constexpr const RungeKuttaNystromCoefficients<N,M> cf;
    constexpr const int Np = N;
    double fs[Np + 1];
    fs[0] = ode(t+cf.c[0]*h, r+v*cf.c[0]*h);
    for (int i=1; i<Np+1; i++) {
      double sum=0e0;
      for (int j=0; j<i; j++) sum += cf.a[i][j]*fs[j];
      fs[i] = ode(t + cf.c[i] * h, r + v * cf.c[i] * h + h*h*sum);
    }
    double sumx = 0e0, sumv = 0e0;
    for (int i=0; i<Np+1; i++) {
      sumx += fs[i] * cf.bhat[i];
      sumv += fs[i] * cf.b[]
    }
    r = r + v*h + 
  }
}// dso

#endif
