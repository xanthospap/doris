#include "egravity.hpp"
#include <cstdio>

int main(int argc, char *argv[]) {
  if (argc != 1) {
    fprintf(stderr, "Usage: %s [GRAVITY MODEL FILE]\n", argv[0]);
    return 1;
  }

  int degree = 4;
  int order = 4;
  double vec[] = {4595216.48810, 2039453.03130, 3912626.78150};
  double Re = 0.6378136460e+07;

#ifdef DEBUG
  Mat2D<MatrixStorageType::RowWise> V(degree + 1, order + 1);
  V.fill_with(0e0);
  Mat2D<MatrixStorageType::RowWise> W(degree + 1, order + 1);
  W.fill_with(0e0);

  int error =
      lagrange_polynomials(vec[0], vec[1], vec[2], Re, degree, order, V, W);
  if (error) {
    fprintf(stderr, "ERROR Computing Lagrange polynomials! #1\n");
    return 1;
  } else {
    printf("Matrix V:\n");
    V.print();
    printf("Matrix W:\n");
    W.print();
  }
#endif

  Mat2D<MatrixStorageType::Trapezoid> Vt(degree + 1, order + 1);
  Mat2D<MatrixStorageType::Trapezoid> Wt(degree + 1, order + 1);
  int error2 =
      lagrange_polynomials(vec[0], vec[1], vec[2], Re, degree, order, Vt, Wt);
  if (error2) {
    fprintf(stderr, "ERROR Computing Lagrange polynomials! #2\n");
    return 1;
  }

#ifdef DEBUG
  printf("Matrix V:\n");
  Vt.print();
  printf("Matrix W:\n");
  Wt.print();
#endif

  int order3 = 2;
  Mat2D<MatrixStorageType::Trapezoid> Vt3(degree + 1, order3 + 1);
  Mat2D<MatrixStorageType::Trapezoid> Wt3(degree + 1, order3 + 1);
  int error3 = lagrange_polynomials(vec[0], vec[1], vec[2], Re, degree, order3,
                                    Vt3, Wt3);
  if (error3) {
    fprintf(stderr, "ERROR Computing Lagrange polynomials! #2\n");
    return 1;
  }

#ifdef DEBUG
  printf("Matrix V:\n");
  Vt3.print();
  printf("Matrix W:\n");
  Wt3.print();
#endif

  return 0;
}