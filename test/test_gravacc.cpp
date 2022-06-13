#include "egravity.hpp"
#include "icgemio.hpp"
#include <cassert>
#include <cstdio>
#include "eigen3/Eigen/Eigen"

using namespace dso;

int main(int argc, char *argv[]) {
  if (argc < 2 || argc > 4) {
    fprintf(stderr,
            "Usage: %s <GRAVITY MODEL FILE> [DEGREE - optional] [ORDER - "
            "optional]\n",
            argv[0]);
    return 1;
  }

  double vec[] = {6525.919e3, 1710.416e3, 2508.886e3};

  // first of all, get the harmonic coefficients and model parameters of from a
  // icgem file
  Icgem gfc(argv[1]);
  // parse the header ...
  if (gfc.parse_header()) {
    fprintf(stderr, "ERROR! Failed to parse icgem header!\n");
    return 1;
  }

  int degree = gfc.degree();
  if (argc >= 3)
    degree = std::atoi(argv[2]);
  printf("Note: Setting degree for spherical harmonics to %d\n", degree);

  int order = degree;
  if (argc == 4)
    order = std::atoi(argv[3]);
  printf("Note: Setting order for spherical harmonics to %d\n",order);
  assert(degree <= gfc.degree() && order <= degree);

  // allocate memory to store harmonic coefficients
  HarmonicCoeffs hc(degree);
  // parse data; store coefficients to hc
  if (gfc.parse_data(degree, order, &hc)) {
    fprintf(stderr, "ERROR! Failed to parse harmonic coefficients\n");
    return 1;
  }

  // de-normalize the harmonics coeffs
  hc.denormalize();

  // ok, we now have the harmonic coefficients; we need the Lagrange
  // polynomials, for the given point
  // Note that we need to compute V and W for degree + 1 (see Montenbruck,
  // 3.2.5)
  // and in case we want partials, we should go up to degree+2 (hence a total
  // of degree+3 elements)
  Mat2D<MatrixStorageType::Trapezoid> V(degree + 3, order + 3),
      W(degree + 3, order + 3);
  if (lagrange_polynomials(vec[0], vec[1], vec[2], gfc.earth_radius(),
                           degree + 2, order + 2, V, W)) {
    fprintf(stderr, "ERROR. Failed to compute Lagrange polynomials\n");
    return 1;
  }

  double acc[3];
  if (grav_potential_accel(degree, order, V, W,
                           hc, acc)) {
    fprintf(stderr, "ERROR. Failed to compute acceleration.\n");
    return 1;
  }

  printf("Acceleration components:\n");
  printf("Xacc: %15.10e\nYacc: %15.10e\nZacc: %15.10e\n", acc[0], acc[1],
         acc[2]);

  // Also compute partials
  Eigen::Matrix<double,3,3> G;
  if (grav_potential_accel(degree, order, V, W,
                           hc, acc, G)) {
    fprintf(stderr, "ERROR. Failed to compute acceleration/partials.\n");
    return 1;
  }
  printf("Acceleration: %+15.9f %+15.9f %+15.9f\n",acc[0], acc[1], acc[2]);
  printf("Partials    : %+15.9f %+15.9f %+15.9f\n",G(0,0),G(0,1),G(0,2));
  printf("            : %+15.9f %+15.9f %+15.9f\n",G(1,0),G(1,1),G(1,2));
  printf("            : %+15.9f %+15.9f %+15.9f\n",G(2,0),G(2,1),G(2,2));

  return 0;
}
