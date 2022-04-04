#ifndef __RUNGE_KUTTA_NYSTROM_ITERGRATOR_HPP__
#define __RUNGE_KUTTA_NYSTROM_ITERGRATOR_HPP__

#include "orbit_integrators.hpp"
#include "orbit_integrators_coefficients.hpp"
#include <cstdio>
#ifdef DEBUG
#include <cassert>
#endif

namespace dso {

template <int N, int M> class RungeKuttaNystromIntegratorCore {
private:
  double h;                        ///< step size
  orbit_integrators::ode2_pv *ode; ///< 2nd degree diff. equation
public:
  constexpr RungeKuttaNystromIntegratorCore(double step_size,
                              orbit_integrators::ode2_pv *function
                              ) noexcept
      : h(step_size), ode(function){};
  
  constexpr double step_size() const noexcept {return step_size;}
  constexpr double& step_size() noexcept {return step_size;}
  
  /// @brief Compute the position and velocity vectors at the next epoch,
  /// aka, t+h
  /// @param[in] t Current epoch
  /// @param[inout] r Solution of position vector at t+h
  /// @param[inout] v Solution of velocity vector at t+h
  /// @return t+h
  double step(double t, Vector3 &r, Vector3 &v) const noexcept {
    constexpr const RungeKuttaNystromCoefficients<N, M> cf;
    Vector3 k[N], dummy;
    ode(t + cf.c[0] * h, r + v * cf.c[0] * h, dummy, k[0]);
    for (int i = 1; i < N; i++) {
      Vector3 sum{{0e0, 0e0, 0e0}};
      for (int j = 0; j < i; j++)
        sum += cf.a[i][j] * k[j];
      ode(t + cf.c[i] * h, r + v * cf.c[i] * h + h * h * sum, dummy, k[i]);
    }

    Vector3 sum{{0e0, 0e0, 0e0}};

    Vector3 rn = r + h * v;
    sum.zero_out();
    for (int i = 0; i < N; i++)
      sum += cf.bhat[i] * k[i];
    rn += (h * h) * sum;

    Vector3 vn = v;
    sum.zero_out();
    for (int i = 0; i < N; i++)
      sum += cf.bdothat[i] * k[i];
    vn += h * sum;

    r = rn;
    v = vn;
    return t + h;
  } 
  /// @brief Compute the position and velocity vectors at the next epoch,
  /// aka, t+h
  /// @param[in] t Current epoch
  /// @param[inout] r Solution of position vector at t+h
  /// @param[inout] v Solution of velocity vector at t+h
  /// @return t+h
  double step_error(double t, Vector3 &r, Vector3 &v, Vector3 &dr) const noexcept {
    constexpr const RungeKuttaNystromCoefficients<N, M> cf;
    Vector3 k[N], dummy;
    ode(t + cf.c[0] * h, r + v * cf.c[0] * h, dummy, k[0]);
    for (int i = 1; i < N; i++) {
      Vector3 sum{{0e0, 0e0, 0e0}};
      for (int j = 0; j < i; j++)
        sum += cf.a[i][j] * k[j];
      ode(t + cf.c[i] * h, r + v * cf.c[i] * h + h * h * sum, dummy, k[i]);
    }

    Vector3 sum{{0e0, 0e0, 0e0}};
    Vector3 drn = r + h*v;
    for (int i=0; i<N; i++) sum += cf.b[i] * k[i];
    drn += h*h*sum;

    Vector3 rn = r + h * v;
    sum.zero_out();
    for (int i = 0; i < N; i++)
      sum += cf.bhat[i] * k[i];
    rn += h * h * sum;

    Vector3 vn = v;
    sum.zero_out();
    for (int i = 0; i < N; i++)
      sum += cf.bdothat[i] * k[i];
    vn += h * sum;

    r = rn;
    v = vn;
    dr = drn;
    return t + h;
  }
};

template <int N, int M, orbit_integrators::StepSizeControl H> class RungeKuttaNystromIntegrator {
private:
  RungeKuttaNystromIntegratorCore<N,M> integrator;

public :
      /// @brief Constructor
      /// @param[in] step_size The step size (aka h)
      /// @param[in] function The 2nd order ODE function; see the
      ///            orbit_integrators::ode2_pv function alias for more info
      /// @note This function ony assigns; to actually initialize the
      /// integrator, you should call the initialize() given initial values.
      constexpr RungeKuttaNystromIntegrator(
          double step_size, orbit_integrators::ode2_pv *function) noexcept
      : integrator(step_size, function) {};
      double step(double t, Vector3 &r, Vector3 &v) const noexcept {return integrator.step(t,r,v); }

}; // RungeKuttaNystromIntegrator

template <int N, int M> class RungeKuttaNystromIntegrator<N,M,orbit_integrators::StepSizeControl::Auto> {
private:
  RungeKuttaNystromIntegratorCore<N,M> integrator;
  double tolerance;                ///< error tolerance for stepsize control
  #ifdef DEBUG
  int max_repetitions = 100;
  int cur_repetition = 0;
  #endif

public:
  /// @brief Constructor
  /// @param[in] step_size The step size (aka h)
  /// @param[in] function The 2nd order ODE function; see the
  ///            orbit_integrators::ode2_pv function alias for more info
  /// @note This function ony assigns; to actually initialize the integrator,
  /// you should call the initialize() given initial values.
  constexpr RungeKuttaNystromIntegrator(double step_size,
                              orbit_integrators::ode2_pv *function,
                              double tlrnc=1e-6) noexcept
      : integrator(step_size, function), tolerance(tlrnc){};

  double step(double t, Vector3 &r, Vector3 &v) noexcept {
    double t0 = t;
    Vector3 r0 = r, v0 = v;
    Vector3 dr;
    double tn = integrator.step_error(t, r, v, dr);
    // error estimate for this step
    double error = (dr - r).norm();
    // new step size ...
    double hnew = 0.9e0 * integrator.h * std::pow(tolerance / error, 1e0 / (N + 1));
    // avoid rapid step-size changes
    if (integrator.step_size()/hnew > 4e0) hnew = 0.25 * integrator.step_size();
    //set new step size
    integrator.step_size() = hnew;
    // repeat step with new step-size if error is larger than tolerance
    if (error>tolerance) {
      return this->step(t0, r0, v0, dr);
    }
    return tn;
  }
}; // RungeKuttaNystromIntegrator

} // namespace dso

#endif
