#include "ode.hpp"
#include <algorithm>

// y <- y + s * a
// where s is a scalar and y, a are vector
double *vecmul(double *y, double s, const double *a, int sz) noexcept {
  for (int i = 0; i < sz; i++)
    y[i] += s * a[i];
  return y;
}

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
int dso::GSOdeSolver::interpolate(double xout, double *yout,
                             double *ypout) noexcept {
  /*
    double x, const double *y, double xout, double *yout,
    double *ypout, int neqn, int kold, const double *phi,
    const double *psi) noexcept {
  */
  if (kold >= 13) {
    fprintf(stderr, "ERROR Invalid call to %s (kold>=13)\n", __func__);
    return 1;
  }

  const double hi = xout - x;
  const int ki = kold + 1;

  // initialize W for computing G
  double w[14];
  for (int i = 1; i < ki + 1; i++)
    w[i] = 1e0 / (double)i;

  // compute G
  double g[14] = {1e0};
  double rho[14] = {1e0};
  double term = 0e0;
  for (int j = 2; j < ki + 1; j++) {
    const double psijm1 = psi(j - 1);
    const double gamma = (hi + term) / psijm1;
    const double eta = hi / psijm1;
    for (int i = 1; i <= ki - j + 1; i++) {
      w[i] = gamma * w[i] - eta * w[i + 1];
    }
    g[j] = w[1];
    rho[j] = gamma * rho[j - 1];
    term = psijm1;
  }

  // set output vectors to zero
  std::fill(yout, yout + neqn, 0e0);
  std::fill(ypout, ypout + neqn, 0e0);

  // interpolate for the solution yout and for the derivative of the
  // solution ypout
  // 𝐲° ← 𝛗 𝐠 (matrix-vector)
  // 𝐲°′ ← 𝛗 𝛒 (matrix-vector)
  for (int j = 1; j < ki + 1; j++) {
    const int i = ki - j + 1;
    const double *PhiColi = phi.slice(i);
    vecmul(yout, g[i], PhiColi, neqn);
    vecmul(ypout, rho[i], PhiColi, neqn);
  }

  // 𝐲° ← 𝐲 + hi 𝐲°
  for (int i = 0; i < neqn; i++)
    yout[i] = yy(i) + hi * yout[i];

  return 0;
}
