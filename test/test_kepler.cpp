#include "astrodynamics.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iersc.hpp"
#include <cstdio>
#include <limits>
#include <random>

constexpr const double pi = iers2010::DPI;
constexpr const double pi2 = iers2010::D2PI;
const double eps_mach = std::numeric_limits<double>::epsilon();
const double eps = 10e0 * eps_mach;
constexpr const int max_it = 20;

std::random_device rd;  // obtain a random number from hardware
std::mt19937 gen(rd()); // seed the generator
std::uniform_real_distribution<> disM(-4e0 * M_PI, 4e0 * M_PI);
// std::uniform_real_distribution<> disM(0e0, 2e0 * M_PI);
std::uniform_real_distribution<> dise(1e-4, 0.999999);

const int num_tests = 100000;

struct MaxDiff {
  double E1{0e0}, E2{0e0}, e, M;
  double operator()() const noexcept { return std::abs(E1 - E2); }
}; // MaxDiff

double Frac(double x) { return x - floor(x); };

double Modulo(double x, double y) { return y * Frac(x / y); }

double EccAnom(double M, double e, int &ok) {

  ok = 0;
  int i = 0;
  double E, f;

  // Starting value
  M = Modulo(M, /*2.0*pi*/ pi2);
  if (e < 0.8)
    E = M;
  else
    E = pi;

  // Iteration
  do {
    f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));
    if (i == max_it) {
      ok = max_it + 1;
      break;
    }
  } while (fabs(f) > eps && ++i < max_it);

#ifdef COUNT_KEPLER_ITERATIONS
  ok = (i >= max_it) ? -1 : i;
#else
  ok = (i >= max_it);
#endif
  return E;
}

unsigned long it1 = 0;
unsigned long it2 = 0;
unsigned long it3 = 0;

int main() {

  MaxDiff md13, md23;
  int ok1, ok2, ok3;
  printf("Note that convergence limit is: %.20e\n", eps);

  // example from Vallado, Example 2-1
  double Mdeg = 235.4e0; // [deg]
  double eval = 0.4e0;   // []
  double VE1 = dso::kepler(eval, dso::deg2rad(Mdeg), ok1, eps);
  double VE2 = dso::kepler_vallado(eval, dso::deg2rad(Mdeg), ok1, eps);
  double VE3 = EccAnom(dso::deg2rad(Mdeg), eval, ok1);
  double RefEdeg = 220.512074767522e0;
  printf("Checking against vallado Example : Ref Result: %.12e [deg]\n",
         RefEdeg);
  printf("                                 : Al1 Result: %.12e Diff: %.12e "
         "[deg]\n",
         dso::rad2deg(VE1), std::abs(dso::rad2deg(VE1) - RefEdeg));
  printf("                                 : Al2 Result: %.12e Diff: %.12e "
         "[deg]\n",
         dso::rad2deg(VE2), std::abs(dso::rad2deg(VE2) - RefEdeg));
  printf("                                 : Al3 Result: %.12e Diff: %.12e "
         "[deg]\n",
         dso::rad2deg(VE3), std::abs(dso::rad2deg(VE3) - RefEdeg));
  assert(std::abs(dso::rad2deg(VE1) - RefEdeg) < 1e-12);
  assert(std::abs(dso::rad2deg(VE2) - RefEdeg) < 1e-12);
  assert(std::abs(dso::rad2deg(VE3) - RefEdeg) < 1e-12);

  for (int i = 0; i < num_tests; i++) {
    double M = disM(gen);
    double e = dise(gen);

    double E1 = dso::kepler(e, M, ok1, eps);
#ifdef COUNT_KEPLER_ITERATIONS
    if (ok1 < 0)
      fprintf(stderr, "#1 Failed to converge!, e=%.5f M=%.5f\n", e, M);
    else
      it1 += ok1;
#else
    if (ok1)
      fprintf(stderr, "#1 Failed to converge!, e=%.5f M=%.5f\n", e, M);
#endif

    double E2 = dso::kepler_vallado(e, M, ok2, eps);
#ifdef COUNT_KEPLER_ITERATIONS
    if (ok2 < 0)
      fprintf(stderr, "#2 Failed to converge!, e=%.5f M=%.5f\n", e, M);
    else
      it2 += ok2;
#else
    if (ok2)
      fprintf(stderr, "#2 Failed to converge!, e=%.5f M=%.5f\n", e, M);
#endif

    [[maybe_unused]] double E3 = EccAnom(M, e, ok3);
#ifdef COUNT_KEPLER_ITERATIONS
    if (ok3 < 0)
      fprintf(stderr, "#3 Failed to converge!, e=%.5f M=%.5f\n", e, M);
    else
      it3 += ok3;
#else
    if (ok3)
      fprintf(stderr, "#3 Failed to converge!, e=%.5f M=%.5f\n", e, M);
#endif

    if (std::abs(E1 - E3) > md13())
      md13 = MaxDiff{E1, E3, e, M};
    if (std::abs(E2 - E3) > md23())
      md23 = MaxDiff{E2, E3, e, M};
  }

  printf("Max difference encountered after %d tests:\n", num_tests);
  printf("Kepler  vs Montenbruck : DeltaEccentricity = %.15e [rad] or %.15e "
         "[deg]\n",
         md13(), dso::rad2deg(md13()));
  printf("Vallado vs Montenbruck : DeltaEccentricity = %.15e [rad] or %.15e "
         "[deg]\n",
         md23(), dso::rad2deg(md23()));
#ifdef COUNT_KEPLER_ITERATIONS
  printf("Average iterations: %lu %lu %lu\n", it1 / num_tests, it2 / num_tests,
         it3 / num_tests);
#endif

  return 0;
}
