#ifndef __DSO_ODE_SOLVERS_HPP__
#define __DSO_ODE_SOLVERS_HPP__

#include "cmat2d.hpp"
#include "enum_bitset.hpp"
#include <cstdint>

namespace dso {

// Function prototype for first order differential equations
// void f (double x, const Vector y, Vector yp[])
typedef void (*ODEfun)(double x,        // Independent variable
                       const double *y, // State vector
                       double *yp,      // Derivative y'=f(x,y)
);

enum class ODEStatusCode : unsigned int {
  INIT = 1,     // Restart integration
  DONE = 2,     // Successful step
  BADACC = 3,   // Accuracy requirement could not be achieved
  NUMSTEPS = 4, // Permitted number of steps exceeded
  STIFF = 5,    // Stiff problem suspected
  INVPARAM = 6  // Invalid input parameters
};              // ODEStatusCode

class GSOdeSolver {
public:
  GSOdeSolver(int num_eqns) noexcept
      : neqn(num_eqns), phi{num_eqns, 17}, VecPool{new double[5 * num_eqns]},
        relerr(0e0), abserr(0e0), kmax(12){};

  ~GSOdeSolver() noexcept {
    if (VecPool)
      delete[] VecPool;
  }
  // C++ implementation
  // https://github.com/xrf/sg-ode/blob/master/sg_ode/ode.c
  // also
  // https://people.math.sc.edu/Burkardt/cpp_src/ode/ode.html
private:
  ODEfun f;                 ///< Differential equation
  int neqn;                 ///< Number of equations
  double MemPool13[5 * 13]; ///< alpha, beta, v, w, psi
  double MemPool14[2 * 14]; ///< sig, g
  dso::Mat2d<dso::MatrixStorageType::ColumnWise> phi; ///<
  double *VecPool;                                    ///< yy, wt, p, yp, ypout
  double h, hold, told, delsgn;
  double relerr; ///< Desired relative accuracy of the solution
  double abserr; ///< Desired absolute accuracy of the solution
  double t;      ///< Value of independent variable
  int ns, k, kold;
  int kmax; // Maximum order
  uint8_t OldPermit, phase1, start, nornd, init;
  bool PermitTOUT; // Flag for integrating past tout
                   // (default = true)

  double &alpha(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 13);
#endif
    return MemPool13[i];
  }
  double &beta(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 13);
#endif
    return MemPool13[13 + i];
  }
  double &v(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 13);
#endif
    return MemPool13[2 * 13 + i];
  }
  double &w(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 13);
#endif
    return MemPool13[3 * 13 + i];
  }
  double &psi(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 13);
#endif
    return MemPool13[4 * 13 + i];
  }
  double &sig(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 14);
#endif
    return MemPool14[i];
  }
  double &g(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 14);
#endif
    return MemPool14[14 + i];
  }
  double &yy(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < neqn);
#endif
    return VecPool[i];
  }
  double &wt(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < neqn);
#endif
    return VecPool[neqn + i];
  }
  double &p(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < neqn);
#endif
    return VecPool[2 * neqn + i];
  }
  double &yp(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < neqn);
#endif
    return VecPool[3 * neqn + i];
  }
  double &ypout(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < neqn);
#endif
    return VecPool[4 * neqn + i];
  }

  static constexpr const double gstr[/*13*/] = {
      0.5e00,    0.0833e0,  0.0417e0,  0.0264e0,  0.0188e0,  0.0143e0, 0.0114e0,
      0.00936e0, 0.00789e0, 0.00679e0, 0.00592e0, 0.00524e0, 0.00468e0};
  static constexpr const double two[/*13*/] = {
      2e0,   4e0,   8e0,    16e0,   32e0,   64e0,  128e0,
      256e0, 512e0, 1024e0, 2048e0, 4096e0, 8192e0};

public:
  // integration step
  int step(double &x, double *y, double &eps, int &crash) noexcept;

  /// interpolation
  int interpolate(double xout, const double *y, double *yout) noexcept;

  EnumBitset<ODEStatusCode> integrate(double &t, double tout, double *y,
                                      int init) noexcept;
};

} // namespace dso

#endif
