#include <algorithm>

/*
INTRP approximates the solution at XOUT by polynomial interpolation.
The methods in STEP approximate the solution near X by a polynomial.
This routine approximates the solution at XOUT by evaluating the
polynomial there.  Information defining this polynomial is passed
from STEP, so INTRP cannot be used alone.

Original FORTRAN77 version by Lawrence Shampine, Marilyn Gordon.

Parameters:

  Input, double X, the point where the solution has been computed.

  Input, double Y[NEQN], the computed solution at X.

  Input, double XOUT, the point at which the solution is desired.

  Output, double YOUT[NEQN], the solution at XOUT.

  Output, double YPOUT[NEQN], the derivative of the solution
  at XOUT.

  Input, int NEQN, the number of equations.

  Input, int KOLD, the order used for the last
  successful step.

  Input, double PHI[NEQN*16], contains information about the
  interpolating polynomial.

  Input, double PSI[12], contains information about the
  interpolating polynomial.

*/
int dso::sg::interp(double x, const double *y, double xout, double *yout,
                     double *ypout, int neqn, int kold, const double *phi,
                     const double *psi) noexcept {
  if (kold >= 13) {
    fprintf(stderr, "ERROR Invalid call to %s (kold>=13)\n", __func__);
    return 1;
  }

  const double hi = xout - x;
  const int ki = kold + 1;

  // initialize W for computing G
  double w[13];
  for (int i=0; i<13; i++) w[i] = 1e0 / (double)(i+1);

  // compute G
  double g[13] = {1e0};
  double rho[13] = {1e0};
  double term = 0e0;
  for (int j = 1; j < ki; j++) {
    const double psijm1 = psi[j - 1];
    const double gamma = (hi + term) / psijm1;
    const double eta = hi / psijm1;
    for (i = 0; i < ki - j; i++) {
      w[i] = gamma * w[i] - eta * w[i+1];
    }
    g[j] = w[0];
    rho[j] = gamma * rho[j - 1];
    term = psijm1;
  }
  
  // set output vectors to zero
  std::fill(yout, yout+neqn, 0e0);
  std::fill(ypout, ypout+neqn, 0e0);

  // interpolate
  for (int i=ki; i>0; i--) {
    const int im1 = i-1;
    for (int k = 0; k < neqn; k++) {
      yout[k] = yout[k] + g[im1] * phi[k + im1 * neqn];
      ypout[k] = ypout[k] + rho[im1] * phi[k + im1 * neqn];
    }
  }

  for (int k=0; k<neqn; k++) yout[k] = y[k] + hi * yout[k];

  return 0;
}
