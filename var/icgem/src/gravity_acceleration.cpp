#include "compact2dmat.hpp"
#include "harmonic_coeffs.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>

void factorial_coeffs(int maxm, double *coeffs) noexcept {
  coeffs[0] = /* 2! / 0! */ 2e0;
  for (int m = 1; m <= maxm; m++)
    coeffs[m] = coeffs[m - 1] * ((m + 2) / m);
  return;
}

/// @ref Montenbruck, Gill, Satellite Orbits, Models Methods Applications;
///      ch. 3.2.5, p. 68
int grav_potential_accel(int order, int degree, double Re, double GM,
                         const Mat2D<MatrixStorageType::Trapezoid> &V,
                         const Mat2D<MatrixStorageType::Trapezoid> &W,
                         const HarmonicCoeffs &hc, double *acc) noexcept {

  double xacc(0e0), yacc(0e0), zacc(0e0);

  // factorial coefficients: (n-m+2)! / (n-m)!
  double *facs = new double[degree + 1];
  factorial_coeffs(degree, facs);

  // m = 0 part
  for (int i = 0; i <= degree; i++) {
    xacc += -hc.C(i, 0) * V(i + 1, 0);
    yacc += -hc.C(i, 0) * W(i + 1, 0);
  }

  double xacc2(0e0), yacc2(0e0);
  // m != 0
  for (int i = 1; i <= degree; i++) {
    for (int j = 1; j <= std::min(i, order); j++) {
      double Cnm = hc.C(i, j);
      double Snm = hc.S(i, j);
      double Vnp1mm1 = V(i + 1, j - 1);
      double Vnp1mp0 = V(i + 1, j);
      double Vnp1mp1 = V(i + 1, j + 1);
      double Wnp1mm1 = W(i + 1, j - 1);
      double Wnp1mp0 = W(i + 1, j);
      double Wnp1mp1 = W(i + 1, j + 1);

      xacc2 += -Cnm * Vnp1mp1 - Snm * Wnp1mp1 +
               facs[i - j] * (Cnm * Vnp1mm1 + Snm * Wnp1mm1);
      yacc2 += -Cnm * Wnp1mp1 + Snm * Vnp1mp1 +
               facs[i - j] * (-Cnm * Vnp1mm1 + Snm * Wnp1mm1);
      zacc += (i - j + 1) * (-Cnm * Vnp1mp0 - Snm * Wnp1mp0);
    }
  }

  delete[] facs;

  xacc += xacc2 / 2e0;
  xacc *= GM / (Re * Re);

  yacc += yacc2 / 2e0;
  yacc *= GM / (Re * Re);

  zacc *= GM / (Re * Re);

  acc[0] = xacc;
  acc[1] = yacc;
  acc[2] = zacc;

  return 0;
}