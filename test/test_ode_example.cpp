#include "eigen3/Eigen/Eigen"
#include "integrators.hpp"

// The ODE system (3 equations)
//   y1'=y2*y3,      y1(0)=0
//   y2'=-y1*y3,     y2(0)=1
//   y3'=-ksq*y1*y2, y3(0)=1
void f([[maybe_unused]] double x,     // Independent variable
       const Eigen::VectorXd &y,      // State (function values)
       Eigen::Ref<Eigen::VectorXd> yp, // Partials/Derivative
       [[maybe_unused]]dso::IntegrationParameters *
       ) noexcept {
  yp(0) = y(1) * y(2);
  // yp(1) = -y(0) * y(2);
  yp(1) = -y(0) * y(2) * 100e0;
  yp(2) = -0.51e0 * y(0) * y(1);
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
  const double relerr_ode = 1e-12;
  const double abserr_ode = 1e-16;

  dso::SGOde Sg(f, neqn, relerr_ode, abserr_ode);

  // initial values
  double t = 0e0;
  Eigen::VectorXd y(3), yy(3);

  // re-initialize
  t = 0e0;
  y << 5e0, 1e0, 1e0;
  printf("- neqn=%3d relerr=%.2e abserr=%.2e iflag=%3d\n", neqn, relerr_ode,
         abserr_ode, 0);
  for (int i = 1; i <= 13; i++) {
    double tout = 5e0 * i;
    Sg.de(t, tout, y, yy);
    printf("# %16.8e %16.8e %16.8e %16.8e\n", t, yy(0), yy(1), yy(2));
    if (Sg.flag() != dso::SGOde::IFLAG::SUCCESS) {
      fprintf(stderr, "Failed to find solution!\n");
      return 1;
    }
    y = yy;
  }
  return 0;
}
