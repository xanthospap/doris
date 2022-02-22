#include "compact2dmat.hpp"
#include "harmonic_coeffs.hpp"
#include "egravity.hpp"
#include <algorithm>
#include <cmath>
#ifdef DEBUG
#include <cstdio>
#include <cassert>
#endif

int dso::grav_potential_accel(
    int degree, int order, double Re, double GM,
    const dso::Mat2D<dso::MatrixStorageType::Trapezoid> &V,
    const dso::Mat2D<dso::MatrixStorageType::Trapezoid> &W,
    const dso::HarmonicCoeffs &hc, double *acc) noexcept {

  double xacc(0e0), yacc(0e0), zacc(0e0);

  // m = 0 part
  for (int i = 0; i <= degree; i++) {
    const double Cn0 = hc.C(i, 0);
    zacc += -(i + 1) * Cn0 * V(i + 1, 0);
    xacc += -Cn0 * V(i + 1, 1);
    yacc += -Cn0 * W(i + 1, 1);
  }

  double xacc2(0e0), yacc2(0e0);
  // m != 0
  for (int i = 1; i <= degree; i++) {
    for (int j = 1; j <= std::min(i, order); j++) {
      const double Cnm = hc.C(i, j);
      const double Snm = hc.S(i, j);
      const double Vnp1mm1 = V(i + 1, j - 1);
      const double Vnp1mp0 = V(i + 1, j);
      const double Vnp1mp1 = V(i + 1, j + 1);
      const double Wnp1mm1 = W(i + 1, j - 1);
      const double Wnp1mp0 = W(i + 1, j);
      const double Wnp1mp1 = W(i + 1, j + 1);
      const double fac = (i - j + 1) * (i - j + 2) / 2e0;

      xacc2 += 0.5e0 * (-Cnm * Vnp1mp1 - Snm * Wnp1mp1) +
               fac * (Cnm * Vnp1mm1 + Snm * Wnp1mm1);
      yacc2 += 0.5e0 * (-Cnm * Wnp1mp1 + Snm * Vnp1mp1) +
               fac * (-Cnm * Vnp1mm1 + Snm * Wnp1mm1);
      zacc += (i - j + 1) * (-Cnm * Vnp1mp0 - Snm * Wnp1mp0);
    }
  }

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
