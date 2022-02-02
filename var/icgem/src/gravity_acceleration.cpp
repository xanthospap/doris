#include "compact2dmat.hpp"
#include "harmonic_coeffs.hpp"
#include <cmath>
#include <cstdio>
#include <algorithm>

    /// @ref Montenbruck, Gill, Satellite Orbits, Models Methods Applications;
    ///      ch. 3.2.5, p. 68
    int grav_potential_accel(int order, int degree,
                             const Mat2D<MatrixStorageType::RowWise> &V,
                             const Mat2D<MatrixStorageType::RowWise> &W,
                             const HarmonicCoeffs &hc) noexcept {

  double xacc(0e0), yacc(0e0), zacc(0e0);

  // m = 0 part
  for (int i = 0; i <= degree; i++) {
    xacc += -hc.C(i, 0) * V(i + 1, 0);
    yacc += -hc.C(i, 0) * W(i + 1, 0);
  }

  // m != 0
  for (int i = 0; i <= degree; i++) {
    for (int j = 1; j <= std::min(i, order); j++) {
        double Cnm = hc.C(i,j);
        double Snm = hc.S(i,j);
        double Vnp1mm1 = V(i+1,j-1);
        double Vnp1mp0 = V(i+1,j);
        double Vnp1mp1 = V(i+1,j+1);
        double Wnp1mm1 = W(i+1,j-1);
        double Wnp1mp0 = W(i+1,j);
        double Wnp1mp1 = W(i+1,j+1);
    }
  }
}