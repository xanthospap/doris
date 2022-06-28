#ifndef __DSO_SGODE_ODE_HPP__
#define __DSO_SGODE_ODE_HPP__

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

class SGOde {
public:
  SGOde(ODEfun _f, int _neqn, double rerr, double aerr, void *_params = nullptr)
      : f(_f), neqn(_neqn), iflag(1), relerr(rerr), abserr(aerr),
        params(_params) {
    Phi  = Eigen::MatrixXd(neqn,16);
    ArraysNeqn = Eigen::MatrixXd(neqn,5); // wt, p, ypout, yp, yy
    Arrays13 = new double[13 * 7];// psi, alpha, beta, v, w, sig, g (col-major)
  }

  ~SGOde() noexcept {
    if (Arrays13) delete[] Arrays13;
  }

  int flag() const noexcept {return iflag;}

  int de(double t0, double tout, const Eigen::VectorXd &y0,
         Eigen::VectorXd &yout) noexcept;
  int step(double &eps, int &crash) noexcept;
  int intrp(double xout, Eigen::VectorXd &yout,
            Eigen::Ref<Eigen::VectorXd> ypout) noexcept;

  Eigen::Ref<Eigen::VectorXd> wt()    noexcept { return ArraysNeqn.col(0); }
  Eigen::Ref<Eigen::VectorXd> p()     noexcept { return ArraysNeqn.col(1); }
  Eigen::Ref<Eigen::VectorXd> yy()    noexcept { return ArraysNeqn.col(2); }
  Eigen::Ref<Eigen::VectorXd> yp()    noexcept { return ArraysNeqn.col(3); }
  Eigen::Ref<Eigen::VectorXd> ypout() noexcept { return ArraysNeqn.col(4); }
  
  double &psi  (int i) noexcept {
  #ifdef DEBUG 
  assert(i>=0 && i<13); 
  #endif 
  return Arrays13[0*13+i]; 
  }
  
  double &alpha(int i) noexcept {
#ifdef DEBUG 
    assert(i >= 0 && i < 13);
#endif 
    return Arrays13[1 * 13 + i];
  }
  
  double &beta(int i) noexcept {
#ifdef DEBUG 
    assert(i >= 0 && i < 13);
#endif 
    return Arrays13[2 * 13 + i];
  }
  
  double *v() noexcept { return Arrays13 + 3 * 13; }
  double &v(int i) noexcept {
#ifdef DEBUG 
    assert(i >= 0 && i < 13);
#endif 
    return Arrays13[3 * 13 + i];
  }
  
  double *w() noexcept { return Arrays13 + 4 * 13; }
  double &w(int i) noexcept {
#ifdef DEBUG 
    assert(i >= 0 && i < 13);
#endif 
    return Arrays13[4 * 13 + i];
  }
  
  double &sig(int i) noexcept {
#ifdef DEBUG 
    assert(i >= 0 && i < 12);
#endif 
    return Arrays13[5 * 13 + i];
  } // rows = 12
  
  double &g(int i) noexcept {
#ifdef DEBUG 
    assert(i >= 0 && i < 12);
#endif 
    return Arrays13[6 * 13 + i];
  } // rows = 12

private:
  ODEfun f;
  int neqn;
  int iflag;
  int start,phase1,nornd,isnold,kold,k;
  Eigen::MatrixXd Phi; // dimension 
  Eigen::MatrixXd ArraysNeqn; // dimension 
  double *Arrays13;
  double h;
  double x;
  double t,told;
  double delsgn;
  double relerr,abserr;
  /// May store a pointer to some king of parameters that are passed in the
  /// ODE function
  void *params{nullptr};
}; // SGOde

} // dso
#endif
