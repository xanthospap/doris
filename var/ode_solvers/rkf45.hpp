#ifndef __DSO_RKF45_ODE_HPP__
#define __DSO_RKF45_ODE_HPP__

#include "eigen3/Eigen/Eigen"

namespace dso {
// Function prototype for first order differential equations
// void f (double x, const Vector y, Vector yp[])
typedef void (*ODEfun)(double x,                 // Independent variable
                       const Eigen::VectorXd &y, // State (function values)
                       /*const Eigen::MatrixXd &Phi,*/
                       Eigen::Ref<Eigen::VectorXd> yp, // Partials/Derivative
                       /*Eigen::MatrixXd &Phip, */
                       void *params);

class RKF45 {
public:
  RKF45(ODEfun _f, int _neqn, double rerr, double aerr, void *_params = nullptr)
      : f(_f), neqn(_neqn), iflag(1), relerr(rerr), abserr(aerr),
        params(_params) {
    F  = Eigen::MatrixXd(neqn,5);
    yp = Eigen::VectorXd(neqn);
  }

private:
  ODEfun f;
  int neqn;
  int iflag, jflag, kflag;
  int nfe,  ///< counter for function evaluations
      kop,  ///< indicator for too many output points
      init; ///< initialization completion indicator
  Eigen::MatrixXd F; // dimension 
  Eigen::VectorXd yp;
  double relerr;
  double abserr;
  double h;
  double savae, savre;
  /// May store a pointer to some king of parameters that are passed in the
  /// ODE function
  void *params{nullptr};

  /// fehlberg fourth-fifth order runge-kutta method
  ///
  /// rkfs integrates a system of first order ordinary differential
  /// equations as described in the comments for rkf45 .
  /// the arrays yp,f1,f2,f3,f4,and f5 (of dimension at least neqn) and
  /// the variables h,savre,savae,nfe,kop,init,jflag,and kflag are used
  /// internally by the code and appear in the call list to eliminate
  /// local retention of variables between calls. accordingly, they
  /// should not be altered. items of possible interest are
  ///       yp - derivative of solution vector at t
  ///       h  - an appropriate stepsize to be used for the next step
  ///       nfe- counter on the number of derivative function evaluations
  int rkfs(double &t, double tout, Eigen::VectorXd &y) noexcept;

  ///  fehlberg fourth-fifth order runge-kutta method
  ///
  /// fehl integrates a system of neqn first order
  /// ordinary differential equations of the form
  ///          dy(i)/dt=f(t,y(1),---,y(neqn))
  /// where the initial values y(i) and the initial derivatives
  /// yp(i) are specified at the starting point t. fehl advances
  /// the solution over the fixed step h and returns
  /// the fifth order (sixth order accurate locally) solution
  /// approximation at t+h in array s(i).
  /// f1,---,f5 are arrays of dimension neqn which are needed
  /// for internal storage.
  /// the formulas have been grouped to control loss of significance.
  /// fehl should be called with an h not smaller than 13 units of
  /// roundoff in t so that the various independent arguments can be
  /// distinguished.
  Eigen::VectorXd fehl(double t, double h, const Eigen::VectorXd &y,
                       const Eigen::VectorXd &yp) noexcept;

public:
  int flag() const noexcept {return iflag;}

  void solve(double &t, double tout, Eigen::VectorXd &y) noexcept {
    if (this->rkfs(t, tout, y) > 0)
      return;
  }

}; // RKF45

} // dso

#endif
