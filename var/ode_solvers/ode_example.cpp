#include "rkf45.hpp"
#include "eigen3/Eigen/Eigen"

// The ODE system (3 equations)
//   y1'=y2*y3,      y1(0)=0
//   y2'=-y1*y3,     y2(0)=1
//   y3'=-ksq*y1*y2, y3(0)=1
void f([[maybe_unused]] double x,                       // Independent variable
       const Eigen::VectorXd &y,       // State (function values)
       Eigen::Ref<Eigen::VectorXd> yp, // Partials/Derivative
       [[maybe_unused]] void *params) {
  yp(0) = y(1) * y(2);
  yp(1) = -y(0) * y(2);
  yp(2) = -0.5e0 * y(0) * y(1);
  return;
}

/*
 *  demonstration program for ode/de/step package.
 * 
 *  reference:shampine and gordon's computer solution of
 *    ordinary differential equations:the initial value problem.
 * 
 *  the package is used to solve the defining equations for the
 *  jacobian elliptic functions:
 *        y1'=y2*y3,      y1(0)=0
 *        y2'=-y1*y3,     y2(0)=1
 *        y3'=-ksq*y1*y2, y3(0)=1
 *  which is solved for ti=i*5 for i=0, 1, ..., 12.
 *  the analytic solution is
 *        y1=sn(t/ksq)
 *        y2=cn(t/ksq)
 *        y3=dn(t/ksq)
 * 
 *  in this case we use ksq=k*k=0.51.
 *  output is written on logical unit lout which has been set to 6.
*/

int main() {
  constexpr const int neqn = 3;
  const double relerr = 1e-9 * 1e2;
  const double abserr = 1e-16 * 1e2;
  const int iflag = 1;

  printf("- neqn=%3d relerr=%.2e abserr=%.2e iflag=%3d\n", neqn, relerr, abserr,
         iflag);
  
  // constructor
  dso::RKF45 rkf45(f, neqn, relerr, abserr);

  // initial values
  double t=0e0;
  Eigen::VectorXd y(3);
  y << 0e0, 1e0, 1e0;
  printf(" %16.8e %16.8e %16.8e %16.8e %3d\n", t, y(0), y(1), y(2), rkf45.flag());

  for (int i=1; i<13; i++) {
    // target time (solution for)
    double tout = 5e0 * i;
    rkf45.solve(t, tout, y);
    printf(" %16.8e %16.8e %16.8e %16.8e %3d\n", t, y(0), y(1), y(2), rkf45.flag());
    // goto (900,11,5,7,9,900),iflag
    //         1  2 3 4 5   6
    switch (rkf45.flag()) {
      case (1):
      case (2):
        ;
        break;
      case (3):
        printf("tolerances too much and have been changed\n");
        break;
      case (4):
        printf("too many steps...eqn. is hard\n");
        return 1;
      case (5):
        printf("eqn. appears to be stiff\n");
        return 1;
      case (6):
        break;
    }
  }
  return 0;
}
