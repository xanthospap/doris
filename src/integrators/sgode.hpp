#ifndef __DSO_SGODE_ODE_HPP__
#define __DSO_SGODE_ODE_HPP__

#include "odefun.hpp"
#include "orbit_integration.hpp"

namespace dso {

class SGOde {
public:
  SGOde(ODEfun _f, int _neqn, double rerr, double aerr, dso::IntegrationParameters *_params = nullptr)
      : f(_f), neqn(_neqn), iflag(1), relerr(rerr), abserr(aerr),
        params(_params) {
    Phi = Eigen::MatrixXd(neqn,16);
    ArraysNeqn = Eigen::MatrixXd(neqn,5); // wt, p, ypout, yp, yy
    Arrays13 = new double[13 * 7]; // psi, alpha, beta, v, w, sig, g (col-major)
  }

  ~SGOde() noexcept {
    if (Arrays13) delete[] Arrays13;
  }

  int flag() const noexcept {return iflag;}
  int &flag() noexcept {return iflag;}

  int de(double &t, double tout, const Eigen::VectorXd &y0,
         Eigen::VectorXd &yout) noexcept;
  int step(double &eps, int &crash) noexcept;
  // ypout is stored in the member variable ypout
  int intrp(double xout, Eigen::VectorXd &yout/*,
            Eigen::Ref<Eigen::VectorXd> ypout*/) noexcept;

  Eigen::Ref<Eigen::VectorXd> wt()    noexcept { return ArraysNeqn.col(0); }
  Eigen::Ref<Eigen::VectorXd> p()     noexcept { return ArraysNeqn.col(1); }
  Eigen::Ref<Eigen::VectorXd> yy()    noexcept { return ArraysNeqn.col(2); }
  // !! Warning !! //
  // Column 3 of ArraysNeqn should be 3, do not change this!
  // see the bug in sgode_step.cpp ~line 330
  Eigen::Ref<Eigen::VectorXd> yp()    noexcept { return ArraysNeqn.col(3); }
  Eigen::Ref<Eigen::VectorXd> ypout() noexcept { return ArraysNeqn.col(4); }
  double &wt   (int i) noexcept { return ArraysNeqn(i,0); }
  double &p    (int i) noexcept { return ArraysNeqn(i,1); }
  double &yy   (int i) noexcept { return ArraysNeqn(i,2); }
  double &yp   (int i) noexcept { return ArraysNeqn(i,3); }
  double &ypout(int i) noexcept { return ArraysNeqn(i,4); }
  
  /// get psi array coefficient (size=12)
  double &psi  (int i) noexcept {
  #ifdef DEBUG 
  assert(i>=0 && i<12); 
  #endif 
  return Arrays13[0*13+i]; 
  }
  
  /// get alpha array coefficient (size=12)
  double &alpha(int i) noexcept {
#ifdef DEBUG 
    assert(i >= 0 && i < 12);
#endif 
    return Arrays13[1 * 13 + i];
  }
  
  /// get beta array coefficient (size=12)
  double &beta(int i) noexcept {
#ifdef DEBUG 
    assert(i >= 0 && i < 12);
#endif 
    return Arrays13[2 * 13 + i];
  }
  
  /// pointer to the first element in the v array (size=12)
  double *v() noexcept { return Arrays13 + 3 * 13; }
  
  /// get v array coefficient (size=12)
  double &v(int i) noexcept {
#ifdef DEBUG 
    assert(i >= 0 && i < 12);
#endif 
    return Arrays13[3 * 13 + i];
  }
  
  /// pointer to the first element in the w array (size=12)
  double *w() noexcept { return Arrays13 + 4 * 13; }
  
  /// get w array coefficient (size=12)
  double &w(int i) noexcept {
#ifdef DEBUG 
    assert(i >= 0 && i < 12);
#endif 
    return Arrays13[4 * 13 + i];
  }
  
  /// get sig array coefficient (size=13)
  double &sig(int i) noexcept {
#ifdef DEBUG 
    assert(i >= 0 && i < 13);
#endif 
    return Arrays13[5 * 13 + i];
  }
  
  /// get g array coefficient (size=13)
  double &g(int i) noexcept {
#ifdef DEBUG 
    assert(i >= 0 && i < 13);
#endif 
    return Arrays13[6 * 13 + i];
  }

public:
//private:
  ODEfun f;
  int neqn;
  int iflag;
  int start,phase1,nornd,isnold,kold,k,ns;
  Eigen::MatrixXd Phi; // dimension 
  Eigen::MatrixXd ArraysNeqn; // dimension 
  double *Arrays13;
  double h,hold;
  double x;
  double /*t,*/told;
  double delsgn;
  double relerr,abserr;
  /// May store a pointer to some king of parameters that are passed in the
  /// ODE function
  dso::IntegrationParameters *params{nullptr};
}; // SGOde

} // dso
#endif
