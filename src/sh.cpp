#include "cmat2d.hpp"
#include "harmonic_coeffs.hpp"
#include <cmath>

/*
 *  Reads:
 * https://boris.unibe.ch/160516/1/adgeo-55-1-2020.pdf
 * https://authors.library.caltech.edu/103344/1/1985SoPh___99__371L.pdf
 * https://space.stackexchange.com/questions/27779/how-can-i-verify-my-reconstructed-gravity-field-of-ceres-from-spherical-harmonic
 * https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/gravity-SphericalHarmonics.pdf
 */

namespace {
inline double _anm(int n, int m) noexcept {
  return std::sqrt(static_cast<double>(4 * n * n - 1) /
                   static_cast<double>(n * n - m * m));
}
}
int legendre_polynomials(
    int maxdegree, double x,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularRowWise> &P) noexcept;

int main(int argc, char *argv[]) {
  const int degree = 120;
  const int order = 120;
  
  dso::HarmonicCoeffs harmonics(degree);
  if (dso::parse_gravity_model(argv[1], degree, order,
                               dso::datetime<dso::nanoseconds>::max(),
                               harmonics, false)) {
    fprintf(stder, "ERROR Failed to parse gravity field\n");
    return 1;
  }

  const Eigen::Matrix<double, 3, 1> r(1.032008992022930295e+06,
                                      -1.537736218408579705e+06,
                                      6.576829252585867420e+06);
}

int sh() noexcept {
}

/*
 * based on "Efficient spherical harmonic transforms aimedat pseudospectral 
 * numerical simulations", Schaeffer N.
 * See also https://arxiv.org/pdf/1410.1748.pdf
 */
int legendre_polynomials(
    int maxdegree, double x,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularRowWise> &P) noexcept {

  double ammf = 1e0;
  const double sqrt1mx2 = std::sqrt(1e0 - x * x);
  const double sqrt14pi = std::sqrt(1e0 / (4e0 * M_PI));
  double pmmfac = 1e0;

  P(0, 0) = 1e0;
  for (int n = 1; n < maxdegree; n++) {
    // Equation (16)
    ammf *= (2e0 * n + 1e0) / (2e0 * n);
    const double ann = sqrt14pi * std::sqrt(ammf);
    // Equation (13)
    P(n, n) = ann * (pmmfac *= sqrt1mx2);
    // Equation (14)
    P(n + 1, n) = _anm(n + 1, n) * P(n, n);
  }

  // missing the last Pnn
  ammf *= (2e0 * maxdegree + 1e0) / (2e0 * maxdegree);
  const double ann = sqrt14pi * std::sqrt(ammf);
  P(maxdegree, maxdegree) = ann * (pmmfac * sqrt1mx2);

  for (int n = 2; n <= maxdegree; n++) {
    const double bn = -std::sqrt(static_cast<double>(2 * n + 1) /
                                 static_cast<double>(2 * n - 3));
    for (int m = 0; m < n - 1; m++) {
      const double bnm =
          bn * std::sqrt(static_cast<double>((n - 1) * (n - 1) - m * m) /
                         static_cast<double>(n * n - m * m));
      P(n, m) = _anm(n, m) * P(n - 1, m) + bnm * P(n - 2, m);
    }
  }

  // all done
  return 0;
}
