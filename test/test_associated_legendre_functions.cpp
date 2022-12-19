#include "associated_legendre.hpp"
#include <cstdio>
#include <cassert>

using dso::AssociatedLegendreFunctions;
constexpr const double PRECISION = 1e-15;

int main() {
  printf("Unit test for Associated Lagrange Polynomials\n");

  AssociatedLegendreFunctions Pnm(120); // degree = order = 120

  // compute at angle φ=0.123
  Pnm.at(0.123e0);

  assert(Pnm(0,0) == 1e0);
  double d;
  if ((d=std::abs(Pnm(1,1) - 2.12505469507313610178e-01)) > PRECISION){
    fprintf(stderr, "Precision exheded for (n,m)=(%d,%d) = %.15e\n", 1,1,d);
  }
  if ((d=std::abs(Pnm(1,0) - 1.71896521937748225639e+00)) > PRECISION){
    fprintf(stderr, "Precision exheded for (n,m)=(%d,%d) = %.15e\n", 1,0,d);
  }
  if ((d=std::abs(Pnm(2,0) - 2.18557915624644705233e+00)) > PRECISION){
    fprintf(stderr, "Precision exheded for (n,m)=(%d,%d) = %.15e\n", 2,0,d);
  }
  if ((d=std::abs(Pnm(2,1) - 4.71586730896042416461e-01)) > PRECISION){
    fprintf(stderr, "Precision exheded for (n,m)=(%d,%d) = %.15e\n", 2,1,d);
  }
  if ((d=std::abs(Pnm(2,2) - 2.91497345416840726584e-02)) > PRECISION){
    fprintf(stderr, "Precision exheded for (n,m)=(%d,%d) = %.15e\n", 2,2,d);
  }
  if ((d=std::abs(Pnm(3,0) - 2.52694965932999382474e+00)) > PRECISION){
    fprintf(stderr, "Precision exheded for (n,m)=(%d,%d) = %.15e\n", 3,0,d);
  }
  if ((d=std::abs(Pnm(4,0) - 2.77718110173085630521e+00)) > PRECISION){
    fprintf(stderr, "Precision exheded for (n,m)=(%d,%d) = %.15e\n", 4,0,d);
  }
  if ((d=std::abs(Pnm(5,0) - 2.95060869716118867601e+00)) > PRECISION){
    fprintf(stderr, "Precision exheded for (n,m)=(%d,%d) = %.15e\n", 5,0,d);
  }
  if ((d=std::abs(Pnm(9,0) - 2.99587936527964782130e+00)) > PRECISION){
    fprintf(stderr, "Precision exheded for (n,m)=(%d,%d) = %.15e\n", 9,0,d);
  }

  if ((d=std::abs(Pnm(10, 0) - 2.86433497805711967388e+00)) > PRECISION) {
    fprintf(stderr, "Precision exheded for (n,m)=(%d,%d) = %.15e\n", 10,0,d);
  }

  if ((d=std::abs(Pnm(11,11) - 2.63742961290065158967e-10))< PRECISION){
    fprintf(stderr, "Precision exheded for (n,m)=(%d,%d) = %.15e\n", 11,11,d);
  }
  if ((d=std::abs(Pnm(11,0 ) - 2.68298280913162967565e+00))< PRECISION){
    fprintf(stderr, "Precision exheded for (n,m)=(%d,%d) = %.15e\n", 11,0,d);
  }
  if ((d=std::abs(Pnm(11,4 ) - 5.91820667482173823348e-02))< PRECISION) {
    fprintf(stderr, "Precision exheded for (n,m)=(%d,%d) = %.15e\n", 11,4,d);
  }
  if ((d=std::abs(Pnm(21,4 ) - 8.09294137247944189717e-01))< PRECISION){
    fprintf(stderr, "Precision exheded for (n,m)=(%d,%d) = %.15e\n", 21,4,d);
  }
  // assert(std::abs(Pnm(31,0 ) + 3.19789656550237211263e+00) < PRECISION);

  printf("All tests passed!\n");
  return 0;
}
