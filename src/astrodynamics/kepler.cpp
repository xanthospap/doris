#include "astrodynamics.hpp"
#include "iers2010/iersc.hpp"
#include <cmath>

constexpr const int max_it = 20;
constexpr const double pi = iers2010::DPI;
constexpr const double pi2 = iers2010::D2PI;

inline double arg02pi(double angle) noexcept {
  const double t = std::fmod(angle, pi2);
  return t + (t<0e0)*pi2;
}

double dso::kepler(double e, double M, int &status,
                   double tolerance_rad) noexcept {
  status = 0;
  M = arg02pi(M);
  const bool elt08 = (e < 8e-1);
  double E = (!elt08) * M_PI + elt08 * M;
  
  double f;
  int it = 0;

  do {
    f = E - e * std::sin(E) - M;
    E = E - f / (1.0 - e * std::cos(E));
  } while (std::abs(f) > tolerance_rad && ++it < max_it);

  status = it >= max_it;
  return E;
}

double dso::kepler_vallado(double e, double M, int &status,
                           double tolerance_rad) noexcept {
  status = 0;
  M = arg02pi(M);
  const bool mr = (M > -M_PI && M < 0) || (M > M_PI);
  const double esign = mr * -e + (!mr) * e;
  double En = M + esign;
  double E;
  int it = 0;

  do {
    E = En;
    En = E + (M - E + e * std::sin(E)) / (1e0 - e * std::cos(E));
  } while (std::abs(E - En) > tolerance_rad && ++it < max_it);

  status = it >= max_it;
  return En;
}

double Frac (double x) { return x-floor(x); };
double Modulo (double x, double y) { return y*Frac(x/y); }
double dso::EccAnom (double M, double e, int &ok)
{

  // Constants
  // const double pi = M_PI;
  const double eps_mach = std::numeric_limits<double>::epsilon();
  const double eps = 100.0*eps_mach;
  ok = 0;

  // Variables

  int    i=0;
  double E, f;

  // Starting value

  M = Modulo(M, /*2.0*pi*/pi2);   
  if (e<0.8) E=M; else E=pi;

  // Iteration
  // printf("\tstarting E=%.15e [M=%+10.4f] (%s)\n", E, M, __func__);

  do {
    f = E - e*sin(E) - M;
    E = E - f / ( 1.0 - e*cos(E) );
    ++i;
    if (i==max_it) {
      // cerr << " convergence problems in EccAnom" << endl;
      ok = 1;
      break;
    }
  }
  while (fabs(f) > eps);

  if (i>=max_it) ok = -1; else ok = i;
  return E;

}
