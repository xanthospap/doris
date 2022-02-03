#include "egravity.hpp"
#include "icgemio.hpp"
#include <cstdio>
#include <cassert>

int main(int argc, char *argv[]) {
  if (argc != 2 && argc != 3) {
    fprintf(stderr, "Usage: %s <GRAVITY MODEL FILE> [DEGREE - optional]\n",
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
  if (argc == 3)
    degree = std::atoi(argv[2]);
  int order = degree;
  assert(degree <= gfc.degree() && order <= degree);

  // allocate memory to store harmonic coefficients
  HarmonicCoeffs hc(degree);
  // parse data; store coefficients to hc
  if (gfc.parse_data(degree, degree, &hc)) {
    fprintf(stderr, "ERROR! Failed to parse harmonic coefficients\n");
    return 1;
  }

  // ok, we now have the harmonic coefficients; we need the Lagrange
  // polynomials.
  // Note that we need to compute V and W for degree + 1 (see Montenbruck,
  // 3.2.5)
  Mat2D<MatrixStorageType::Trapezoid> V(degree + 2, order + 2),
      W(degree + 2, order + 2);
  if (lagrange_polynomials(vec[0], vec[1], vec[2], gfc.earth_radius(), degree + 1,
                           order + 1, V, W)) {
    fprintf(stderr, "ERROR. Failed to compute Lagrange polynomials\n");
    return 1;
  }
  #ifdef DEBUG
  //for (int i=0; i<=degree+1; i++) {
  //    for (int j=0; j<=i; j++) {
  //        printf("%+15.10e ", V(i,j));
  //    }
  //    printf("\n");
  //}
  //for (int i=0; i<V.num_elements(); i++) printf("%+15.10e ", *(V.data()+i));
  //printf("\n");
  #endif
  double acc[3];
  if (grav_potential_accel(degree, order, gfc.earth_radius(), gfc.gm(), V, W, hc, acc)) {
      fprintf(stderr, "ERROR. Failed to compute acceleration.\n");
      return 1;
  }

  printf("Acceleration components:\n");
  printf("Xacc: %15.10e\nYacc: %15.10e\nZacc: %15.10e\n", acc[0], acc[1], acc[2]);

  return 0;
}