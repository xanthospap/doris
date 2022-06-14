#include "cmat2d.hpp"
#include "egravity.hpp"
#include "factorial.hpp"
#include "harmonic_coeffs.hpp"
#include <algorithm>
#include <cmath>
#ifdef DEBUG
#include <cassert>
#include <cstdio>
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

int dso::grav_potential_accel(
    int degree, int order, double Re, double GM,
    const dso::Mat2D<dso::MatrixStorageType::Trapezoid> &V,
    const dso::Mat2D<dso::MatrixStorageType::Trapezoid> &W,
    const dso::HarmonicCoeffs &hc, double *acc,
    Eigen::Matrix<double, 3, 3> &partials) noexcept {

  // factorial look-up table up to N=120
  // constexpr const dso::FactorialLookUpTable<long unsigned, 120> FTable;
  // static_assert((double)FTable.factorial(1) == 1e0);

  // acceleration
  double xacc(0e0), yacc(0e0), zacc(0e0);
  // partials
  double daxdx_m0(0e0), daxdy_m0(0e0), daxdz_m0(0e0), daydz_m0(0e0),
      dazdz_m0(0e0);

  double fac = 2e0; /* recursion of (n+2)! / n! */
  // m = 0 part
  for (int i = 0; i <= degree; i++) {
    const double Cn0 = hc.C(i, 0);
    /*const double Sn0 = hc.S(i, 0) = 0e0;*/

    // acceleration
    zacc += -(i + 1) * Cn0 * V(i + 1, 0);
    xacc += -Cn0 * V(i + 1, 1);
    yacc += -Cn0 * W(i + 1, 1);

    // partials ...
    daxdx_m0 += 0.5e0 * (Cn0 * V(i + 2, 2) - fac * (Cn0 * V(i + 2, 0)));
    daxdy_m0 += 0.5e0 * (Cn0 * W(i + 2, 2));
    daxdz_m0 += (i + 1) * Cn0 * V(i + 2, 1);
    daydz_m0 += (i + 1) * Cn0 * W(i + 2, 1);
    dazdz_m0 += (double)((i + 1) * (i + 2)) *
                (Cn0 * V(i + 2, 0) /*+ Sn0 * W(i + 2, 0)*/);

    fac *= ((double)(i + 3) / (double)(i + 1));
  }

  double daxdx_m1(0e0), daxdy_m1(0e0), dazdz_m1(0e0);
  // m = 1 part (only considered for partials)
  fac = 2e0; /* recursion of (n+1)! / (n-1)! */
  for (int i = 1; i <= degree; i++) {
    const double Cn1 = hc.C(i, 1);
    const double Sn1 = hc.S(i, 1);
    daxdx_m1 += 0.25e0 * (Cn1 * V(i + 2, 3) + Sn1 * W(i + 2, 3) +
                          fac * (-3e0 * Cn1 * V(i + 2, 1) - Sn1 * W(i + 2, 1)));
    daxdy_m1 += 0.25e0 * (Cn1 * W(i + 2, 3) - Sn1 * V(i + 2, 3) +
                          fac * (-Cn1 * W(i + 2, 1) - Sn1 * V(i + 2, 1)));
    dazdz_m1 += (double)(i * (i + 1)) * (Cn1 * V(i + 2, 1) + Sn1 * W(i + 2, 1));

    fac *= ((double)(i + 2) / (double)i);
  }

  double xacc2(0e0), yacc2(0e0);
  double daxdx_m2(0e0), daxdy_m2(0e0), daxdz_m2(0e0), daydz_m2(0e0),
      dazdz_m2(0e0);
  // m > 1
  for (int i = 1; i <= degree; i++) {
    double fac2, fac3, fac4;
    if (i > 1) {
      /* (n-m+2)! / (n-m)! for n=2,... , starting with n=i,m=2 */
      fac2 = dso::factorialRatio<double>(i - 2 + 2, i - 2);
      //(double)FTable.factorial(i - 2 + 2) / (double)FTable.factorial(i - 2);
      /* (n-m+3)! / (n-m)! for n=2,... */
      fac3 = dso::factorialRatio<double>(i - 2 + 3, i - 2);
      //(double)FTable.factorial(i - 2 + 3) / (double)FTable.factorial(i - 2);
      /* (n-m+4)! / (n-m)! for n=2,... */
      fac4 = dso::factorialRatio<double>(i - 2 + 4, i - 2);
      //(double)FTable.factorial(i - 2 + 4) / (double)FTable.factorial(i - 2);
    }

    for (int j = 1; j <= std::min(i, order); j++) {
      const double Cnm = hc.C(i, j);
      const double Snm = hc.S(i, j);
      const double Vnp1mm1 = V(i + 1, j - 1);
      const double Vnp1mp0 = V(i + 1, j);
      const double Vnp1mp1 = V(i + 1, j + 1);
      const double Wnp1mm1 = W(i + 1, j - 1);
      const double Wnp1mp0 = W(i + 1, j);
      const double Wnp1mp1 = W(i + 1, j + 1);
      const double afac = (i - j + 1) * (i - j + 2) / 2e0;

      // acceleration
      xacc2 += 0.5e0 * (-Cnm * Vnp1mp1 - Snm * Wnp1mp1) +
               afac * (Cnm * Vnp1mm1 + Snm * Wnp1mm1);
      yacc2 += 0.5e0 * (-Cnm * Wnp1mp1 + Snm * Vnp1mp1) +
               afac * (-Cnm * Vnp1mm1 + Snm * Wnp1mm1);
      zacc += (i - j + 1) * (-Cnm * Vnp1mp0 - Snm * Wnp1mp0);

      // partials ... (already computed for m==1)
      if (j > 1) {
        daxdx_m2 +=
            0.25e0 * (Cnm * V(i + 2, j + 2) + Snm * W(i + 2, j + 2) +
                      2e0 * fac2 * (-Cnm * V(i + 2, j) - Snm * W(i + 2, j)) +
                      fac4 * (Cnm * V(i + 2, j - 2) + Snm * W(i + 2, j - 2)));
        daxdy_m2 +=
            0.25e0 * (Cnm * W(i + 2, j + 2) - Snm * V(i + 2, j + 2) +
                      fac4 * (-Cnm * W(i + 2, j - 2) + Snm * V(i + 2, j - 2)));

        daxdz_m2 +=
            (static_cast<double>(i - j + 1) / 2e0) *
                (Cnm * V(i + 2, j + 1) + Snm * W(i + 2, j + 1)) +
            fac3 * 0.5e0 * (-Cnm * V(i + 2, j - 1) - Snm * W(i + 2, j - 1));

        daydz_m2 +=
            (static_cast<double>(i - j + 1) / 2e0) *
                (Cnm * W(i + 2, j + 1) - Snm * V(i + 2, j + 1)) +
            fac3 * 0.5e0 * (Cnm * W(i + 2, j - 1) - Snm * V(i + 2, j - 1));

        dazdz_m2 += fac2 * (Cnm * V(i + 2, j) + Snm * W(i + 2, j));

        // printf("\tUsed factorial values: %.2f %.2f %.2f for n=%d and m=%d\n",
        // fac2, fac3, fac4, i, j);
        fac2 *= ((double)(i - j) / (double)(i - j + 2));
        fac3 *= ((double)(i - j) / (double)(i - j + 3));
        fac4 *= ((double)(i - j) / (double)(i - j + 4));
      }
    }
  }

  // Sum-up acceleration
  xacc += xacc2;
  xacc *= GM / (Re * Re);

  yacc += yacc2;
  yacc *= GM / (Re * Re);

  zacc *= GM / (Re * Re);

  acc[0] = xacc;
  acc[1] = yacc;
  acc[2] = zacc;

  // Sum-up partial derivatives
  // printf("\tdx/dx = %.9f + %.9f + %.9f\n", daxdx_m0, daxdx_m1, daxdx_m2);
  // printf("\tdx/dy = %.9f + %.9f + %.9f\n", daxdy_m0, daxdy_m1, daxdy_m2);
  // printf("\tdx/dz = %.9f + %.9f \n", daxdz_m0, daxdz_m2);
  partials(0, 0) = daxdx_m0 + daxdx_m1 + daxdx_m2; // dax/dx
  partials(0, 1) = daxdy_m0 + daxdy_m1 + daxdy_m2; // dax/dy
  partials(0, 2) = daxdz_m0 + daxdz_m2;            // dax/dz

  partials(1, 0) = partials(0, 1);      // day/dx
  partials(1, 1) = 1e0;                 // day/dy
  partials(1, 2) = daydz_m0 + daydz_m2; // day/dz

  partials(2, 0) = partials(0, 2);                 // daz/dx
  partials(2, 1) = partials(1, 2);                 // daz/dy
  partials(2, 2) = dazdz_m0 + dazdz_m1 + dazdz_m2; // daz/dz

  // computation of day / dy from the fact that:
  // dax / dx + day / dy + daz / dz = 0
  partials(1, 1) = -(partials(0, 0) + partials(2, 2));

  partials *= (GM / Re / Re / Re);

  return 0;
}
