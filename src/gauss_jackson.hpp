#ifndef __GAUSS_JACKSON_ORBIT_ITERGRATOR_HPP__
#define __GAUSS_JACKSON_ORBIT_ITERGRATOR_HPP__

#include "iers2010/matvec.hpp"
#include <type_traits>
#include <cstdio>
#ifdef DEBUG
#include <cassert>
#endif

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

/// @brief Valid order for integrators; check at compile-time
#if __cplusplus >= 202002L
    template<int N> concept ValidOrbitIntegratorOrder = N > 1 && N < 9;
#else
    template<int N> struct ValidOrbitIntegratorOrder {
        static bool const value = N > 1 && N < 9;
    };
#endif
} // namespace orbit_integrators

/// @brief Integrate a 2nd order ODE via Runge-Kutta 4
/// Compute the solution of a 2nd order ODE using the Runge-Kutta 4 method. The 
/// 2nd order ODE is passed in as a function pointer; the input arguments (in
/// corresponding order) are: the time t, the position vector (at t), the 
/// velocity vector (at time t) and the output acceleration vector, aka the
/// ODE has a signature of:
/// void f(double t, const Vector3 &pos, const Vector3 &vel, Vector3 &acc)
/// noexcept
/// The function will compute the (approximate) solution (position and velocity) 
/// at t+h
/// @param[in] t Starting time
/// @param[in] h RK4 step size
/// @param[inout] r At input, position vector at t; at output RK4 solution, aka
///               position vector (solution) at t+h
/// @param[inout] v At input, velocity vector at t; at output RK4 solution, aka
///               velocity vector (solution) at t+h
/// @param[in] ode Pointer to the 2nd order ODE function
/// @return t+h Aka the time at the solution epoch
double rk4_integrate(double t, double h, Vector3 &r, Vector3 &v,
                     orbit_integrators::ode2_pv *ode) noexcept {
  Vector3 ka1, ka2, ka3, ka4;

  const Vector3 r1 = r;
  const Vector3 v1 = v;
  ode(t, r1, v1, ka1);

  const Vector3 r2 = r + v1 * (h / 2e0);
  const Vector3 v2 = v + ka1 * (h / 2e0);
  ode(t, r2, v2, ka2);

  const Vector3 r3 = r + v2 * (h / 2e0);
  const Vector3 v3 = v + ka2 * (h / 2e0);
  ode(t, r3, v3, ka3);

  const Vector3 r4 = r + v3 * h;
  const Vector3 v4 = v + ka3 * h;
  ode(t, r4, v4, ka4);

  r += (v1 + 2e0 * v2 + 2e0 * v3 + v4) * (h / 6e0);
  v += (ka1 + 2e0 * ka2 + 2e0 * ka3 + ka4) * (h / 6e0);
  return t + h;
}

#if __cplusplus >= 202002L
template <int Order>
class GaussJacksonIntegrator requires orbit_integrators::
    ValidOrbitIntegratorOrder {
#else
template <int Order,
          typename = std::enable_if_t<
              orbit_integrators::ValidOrbitIntegratorOrder<Order>::value>>
class GaussJacksonIntegrator {
#endif
private:
  double h;                        ///< step size
  orbit_integrators::ode2_pv *ode; ///< 2nd degree diff. equation
  Vector3 D[Order + 1]; ///< Backward differences of acceleration at t
  Vector3 d[Order + 1]; ///< Backward differences of acceleration at t+h
  Vector3 S1;           ///< First sum of acceleration
  Vector3 S2;           ///< Second sum of acceleration

public:
  /// @brief Constructor
  /// @param[in] step_size The step size (aka h)
  /// @param[in] function The 2nd order ODE function; see the 
  ///            orbit_integrators::ode2_pv function alias for more info
  /// @note This function ony assigns; to actually initialize the integrator, 
  /// you should call the initialize() given initial values.
  GaussJacksonIntegrator(double step_size, 
                         orbit_integrators::ode2_pv *function) noexcept
      : h(step_size), ode(function){};

  /// @brief Compute the position and velocity vectors at the next epoch,
  /// aka, t+h
  /// @param[in] t Current epoch
  /// @param[in] r Solution of position vector at t+h
  /// @param[in] r Solution of velocity vector at t+h 
  /// @return t+h
  double step(double t, Vector3 &r, Vector3 &v) noexcept {
    constexpr const double gp[] = {1e0,
                                   0.5e0,
                                   5e0 / 12e0,
                                   3e0 / 8e0,
                                   251e0 / 720e0,
                                   95e0 / 288e0,
                                   19087e0 / 60480e0,
                                   5257e0 / 17280e0,
                                   1070017e0 / 3628800e0};
    constexpr const double dp[] = {1e0,
                                   0e0,
                                   1e0 / 12e0,
                                   1e0 / 12e0,
                                   19e0 / 240e0,
                                   3e0 / 40e0,
                                   863e0 / 12096e0,
                                   275e0 / 4032e0,
                                   33953e0 / 518400e0};

    // predictor
    Vector3 rp = dp[0] * S2;
    for (int i = 2; i < Order + 2; i++)
      rp += dp[i] * D[i - 2];
    rp *= (h * h);
    Vector3 vp = gp[0] * S1;
    for (int i = 1; i < Order + 1; i++)
      vp += gp[i] * D[i - 1];
    vp *= h;

    // Update backwards difference table
    ode(t + h, rp, vp, d[0]);
    for (int i = 1; i < Order; i++)
      d[i] = d[i - 1] - D[i - 1];
    for (int i = 0; i < Order; i++)
      D[i] = d[i];
    S1 += d[0];
    S2 += S1;
    
    // Update independent variable and solution
    r = rp;
    v = vp;
    return t+h;
  }
  
  /// @brief Initialize the integrator, using RK4.
  /// Compute the acceleration at Order-1 past times (including at time t0).
  /// Use them to compute the backwards differences
  /// @param[in] t0 Starting epoch; the function will compute approximate
  ///            accelerations at epochs t-i*h, i=1,...,Order-1 using the
  ///            solution obtained by an RK4
  /// @param[in] r0 Intial position vector (at t=t0)
  /// @param[in] vo Initial velocity vector (at t=t0)
  void initialize(double t0, const Vector3 &r0, const Vector3 &v0) noexcept {
    // Create table of accelerations at past times t-3h, t-2h, and t-h using
    // RK4 steps
    Vector3 rt(r0), vt(v0);
    double t = t0;

    // compute accelaration at t0
    ode(t0, r0, v0, D[0]);

    // use RK4 to get an approximate solution for t-i*h, i=1,...,Order-1
    // use the solution to compute the acceleration at the respective epochs
    for (int i = 1; i < Order; i++) {
      t = rk4_integrate(t, -h, rt, vt, ode);
      ode(t, rt, vt, D[i]);
    }

    // Compute backwards differences
    // assuming an order of 4, this will yield (D is the backward difference
    // operator):
    // D[0] = a_t0
    // D[1] = D a_t0
    // D[2] = D^2 a_t0
    // D[3] = D^3 a_t0
    for (int order = 1; order < Order; order++)
      for (int j = Order - 1; j >= order; j--)
        D[j] = D[j - 1] - D[j];

    constexpr const double SCdeltaCoeffs[] = {1e0,
                                              -1e0,
                                              1e0 / 12e0,
                                              0e0,
                                              -1e0 / 240e0,
                                              -1e0 / 240e0,
                                              -221e0 / 60480e0,
                                              -19e0 / 60480e0,
                                              -9829e0 / 3628800e0};

    constexpr const double AMgammaCoeffs[] = {1e0,
                                              -0.5e0,
                                              -1e0 / 12e0,
                                              -1e0 / 24e0,
                                              -19e0 / 720e0,
                                              -3e0 / 160e0,
                                              -863e0 / 60480e0,
                                              -275e0 / 24192e0,
                                              -33953e0 / 3628800e0};

    // Initialize backwards sums using GJ corrector via Adams—Moulton and 
    // Cowell formulas (see Monentbruck, Eq. 4.99)
    S1 = v0 / h;
    for (int i = 1; i <= Order; i++)
      S1 -= AMgammaCoeffs[i] * D[i - 1];
    S2 = r0 / (h * h) - SCdeltaCoeffs[1] * S1;
    for (int i = 2; i <= Order + 1; i++)
      S2 -= SCdeltaCoeffs[i] * D[i - 2];

    // note S1 = D^(-1) a_t0
    //      S2 = D^(-2) a_t0
    /*printf("\tS1=[%.5f, %.5f, %.5f]\n", S1.data[0], S1.data[1], S1.data[2]);
    printf("\tS2=[%.5f, %.5f, %.5f]\n", S2.data[0], S2.data[1], S2.data[2]); */
#ifdef DEBUG
    Vector3 nullvec{{0e0,0e0,0e0}};
    for (int i=0; i<Order+1; i++) assert(d[i] == nullvec);
#endif
  }
}; // GaussJacksonIntegrator

} // namespace dso

#endif
