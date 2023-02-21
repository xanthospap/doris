#ifndef __DSO_DDEABM_ODE_HPP__
#define __DSO_DDEABM_ODE_HPP__

#include "odefun.hpp"
#include "orbit_integration.hpp"

/* see https://github.com/jacobwilliams/ddeabm/blob/master/src/ddeabm_module.F90 */

namespace dso {

class Ddeabm {
  /// The constant maxnum is the maximum number of steps allowed
  static constexpr const int maxnum = 500;
  static constexpr const int maxorder = 12;
private:
  /// slope function (uses params pointer for evaluation)
  ODEfun f; 
  ///  number of equations (aka practically size of arrays used)
  int neqn; 
  double relerr, abserr; // rtol atol
  double trelerr, tabserr; // rtol_tmp, atol_tmp
  bool scalar_tols {true};
  int info[4] = {0,0,0,0};
  int initial_step_mode {1};
  double initial_step_size = 0e0;
  int icount{0};
  double tprev=0e0;
  double h;
  double eps;
  double x;
  double xold;
  double hold;
  double told;
  double delsgn;
  double tstop;
  int start, phase1, nornd,stiff,intout;
  int ns,kord,init,ksteps,kle4,iquit,kprev,ivc,iv[10],kgi,error;
  Eigen::VectorXd ypout, yp, yy, wt, p;
  Eigen::Matrix<double,13,8> // alpha, beta, psi, v, w, sig, g, gi
  Eigen::MatrixXd Phi;
  double gstr[maxorder+2];

  // Compute the array GSTR as in page 159 of Shampine & Gordon, 1975. These 
  // are the γ^{*}_{i} coefficients, so that 
  //           GSTR(I) = | γ^{*}_{I} | for i=[0,1,2,...,MAXORDER+1]
  // Obviously the zero index is never used in FORTRAN.
  // Originally, this function computed powers of two (array named TWO), but 
  // this is very cheap to care (just use ints and cast them).
  void dsteps() {
    gstr[0] = 1e0;
    for (int j=1; j<maxorder+2; j++) {
      gstr[j] = 0e0;
      for (int i=0; i<j; i++) {
        gstr[j] -= gstr[i] / (j+1-i);
      }
      gstr[j] = -gstr[j];
    }
    return;
  }// dsteps
};

}//namespace dso

#endif
