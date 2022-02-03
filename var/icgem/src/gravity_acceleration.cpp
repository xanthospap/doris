#include "compact2dmat.hpp"
#include "harmonic_coeffs.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cassert>

double normfac(int n, int m) noexcept {
  static unsigned long facar[] = {1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,6227020800,87178291200,1307674368000,20922789888000,355687428096000,6402373705728000,121645100408832000};
  assert(n+m < (int)sizeof(facar) && n-m>=0);
  double f = (double)facar[n+m] / (double)facar[n-m];
  return std::sqrt(f/(2-(m==0))/(2*n+1));
}

/// @ref Montenbruck, Gill, Satellite Orbits, Models Methods Applications;
///      ch. 3.2.5, p. 68
int grav_potential_accel(int degree, int order, double Re, double GM,
                         const Mat2D<MatrixStorageType::Trapezoid> &V,
                         const Mat2D<MatrixStorageType::Trapezoid> &W,
                         const HarmonicCoeffs &hc, double *acc) noexcept {

  double xacc(0e0), yacc(0e0), zacc(0e0);

  printf("ax=");
  // m = 0 part
  for (int i = 0; i <= degree; i++) {
    const double nfac = normfac(i,0);
    const double Cn0 = hc.C(i,0)/nfac;
    zacc += (i+1) * Cn0 * V(i+1,0);
    xacc += -Cn0 * V(i + 1, 1);
    yacc += -Cn0 * W(i + 1, 1);
    printf("-%15.10e * %15.10e\n", Cn0, V(i + 1, 1));
  }

  double xacc2(0e0), yacc2(0e0);
  // m != 0
  for (int i = 1; i <= degree; i++) {
    for (int j = 1; j <= std::min(i, order); j++) {
      const double nfac = normfac(i,j);
      const double Cnm = hc.C(i, j) / nfac;
      const double Snm = hc.S(i, j) / nfac;
      const double Vnp1mm1 = V(i + 1, j - 1);
      const double Vnp1mp0 = V(i + 1, j);
      const double Vnp1mp1 = V(i + 1, j + 1);
      const double Wnp1mm1 = W(i + 1, j - 1);
      const double Wnp1mp0 = W(i + 1, j);
      const double Wnp1mp1 = W(i + 1, j + 1);
      const double fac = (i-j+1) * (i-j+2) / 2e0;

      xacc2 += 0.5e0 * (-Cnm * Vnp1mp1 - Snm * Wnp1mp1) +
               fac * (Cnm * Vnp1mm1 + Snm * Wnp1mm1);
      printf("+0.5 * (-%15.10e * %15.10e - %15.10e * %15.10e)\n", Cnm, Vnp1mp1, Snm, Wnp1mp1);
      printf("+%15.10e * (%15.10e * %15.10e + %15.10e * %15.10e)\n", fac, Cnm, Vnp1mm1, Snm, Wnp1mm1);
      yacc2 += 0.5e0 * (-Cnm * Wnp1mp1 + Snm * Vnp1mp1) +
               fac * (-Cnm * Vnp1mm1 + Snm * Wnp1mm1);
      zacc += (i - j + 1) * (-Cnm * Vnp1mp0 - Snm * Wnp1mp0);
    }
  }
  printf("\n");

  xacc += xacc2;
  xacc *= GM / (Re * Re);

  yacc += yacc2;
  yacc *= GM / (Re * Re);

  zacc *= GM / (Re * Re);

  acc[0] = xacc;
  acc[1] = yacc;
  acc[2] = zacc;

  return 0;
}