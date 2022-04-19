#include "astrodynamics.hpp"
#include "geodesy/units.hpp"
#include <limits>
#include <random>
#include <cstdio>

std::random_device rd;  // obtain a random number from hardware
std::mt19937 gen(rd()); // seed the generator
//std::uniform_real_distribution<> disM(-4e0 * M_PI, 4e0 * M_PI);
std::uniform_real_distribution<> disM(0e0, 2e0 * M_PI);
std::uniform_real_distribution<> dise(1e-4, 0.999999);

const int num_tests = 1000;

struct MaxDiff {
  double E1{0e0},E2{0e0},e,M;
  double operator()() const noexcept { return std::abs(E1-E2); }
}; // MaxDiff

unsigned long it1 = 0;
unsigned long it2 = 0;
unsigned long it3 = 0;

int main() {

  MaxDiff md13, md23;
  int ok1, ok2, ok3;
  for (int i=0; i<num_tests; i++) {
    double M = disM(gen);
    double e = dise(gen);
    
    //printf("Testing with e=%.5f M=%.5f, M1=%+10.5f, M2=%+10.5f\n", e,M, Modulo(M, 2.0*M_PI), n(M));
    double E1 = dso::kepler(e, M, ok1);
    if (ok1<0)
      fprintf(stderr, "#1 Failed to converge!, e=%.5f M=%.5f\n", e, M);
    else
      it1 += ok1;
    double E2 = dso::kepler_vallado(e, M, ok2);
    if (ok2<0)
      fprintf(stderr, "#2 Failed to converge!, e=%.5f M=%.5f\n", e, M);
    else
      it2 += ok2;
    [[maybe_unused]]double E3 = dso::EccAnom(M, e, ok3);
    if (ok3<0)
      fprintf(stderr, "#3 Failed to converge!, e=%.5f M=%.5f\n", e, M);
    else
      it3 += ok3;
    printf("Iterations %d %d %d\n", ok1, ok2, ok3);
    
    if (std::abs(E1-E3) > md13())
      md13 = MaxDiff{E1,E3,e,M};
    if (std::abs(E2-E3) > md23())
      md23 = MaxDiff{E2,E3,e,M};
  }

  printf("Max difference encountered after %d tests:\n", num_tests);
  printf("Kepler  vs Montenbruck : DeltaEccentricity = %.15e [rad] or %.15e [deg]\n", md13(), dso::rad2deg(md13()));
  printf("Vallado vs Montenbruck : DeltaEccentricity = %.15e [rad] or %.15e [deg]\n", md23(), dso::rad2deg(md23()));
  printf("Average iterations: %lu %lu %lu\n", it1/num_tests, it2/num_tests, it3/num_tests);

  return 0;
}
