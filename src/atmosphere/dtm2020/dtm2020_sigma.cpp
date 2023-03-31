#include "atmosphere/dtm2020/dtm2020.hpp"

/* This is actually not used at the time, as it provides an one-sigma error 
 * evaluation module which is not needed for now.
 * Avoid compilation.
 */
#ifdef USE_DTM2020_SIGMA_EVALUATION
#include <array>

/*
 * This is a translated version of the FORTAN source code found at
 * https://github.com/swami-h2020-eu/mcm/blob/main/src/dtm2020/sigma_function.f90
 *
 * Original (FORTRAN) Author and translated Version:
 * Author : Claude Boniface (CNES)
 * last version updated : 05/01/2021
 */

namespace {
/* general polynomial function of compile time order
 * usage:
 *   Poly<3> p = {1e0, 2e0, 3e0}; p(5);
 *   double p1 = Poly<3>({1,2,3})(5); // = 1*5^0 + 2*5^1 + 3*5^2
 */
template <int N> struct Poly {
  /* c0, c1, c2, ..., cn */
  std::array<double, N> mc;
  /* c0 + c1*x + c2*x^2 + ... + cn*x^n */
  constexpr double operator()(double x) noexcept {
    double s = 0e0;
    double t = 1e0;
    for (int i = 0; i < N; i++) {
      s += mc[i] * t;
      t *= x;
    }
    return s;
  }
}; // struct Poly

/* matching fortran function std_FP_M(x) */
double solar_flux(double x) noexcept {
  return Poly<2>({26.04e0, -3.36e-2})(x);
}
/* matching fortran function std_KP_M(x) */
double magnetic_index(double x) noexcept {
  return Poly<4>({22.47e0, 7.46e-1, -1.99e-1, 6.56e-2})(x);
}
/* matching function std_alat_M(x) */
double latitude(double x) noexcept {
  return Poly<3>({22.33e0,-7.19e-3,4.72e-4})(x);
}
/* matching std_hl_M */
double local_hour(double x) noexcept {
  return Poly<4>({28.53e0,-3.95e-1,-5.64e-2,3.15e-3})(x);
}
/* matching std_alti_M */
double altitude(double x) noexcept {
  return Poly<2>({-9.73e0,7.92e-2})(x);
}
/* matching function std_day_M */
double day_year(double x) noexcept {
  return Poly<2>({23.15e0,-3.91e-3})(x);
}

}// unnamed namespace

int dso::Dtm2020::sigma_function() noexcept {
  if (altitude() < 200e0) {
    /* Global 1-sigma calculation by means of statistical binning of data
     * (GOCE, CHAMP, GRACE and SWARM) and Hyperplan equation (6 dimensions)
     * Polynomial extrapolation by fixing altitude value at 200 km */
    Copt_L = 0.929e0;
    Res_L = 112.79e0;
    sigma_alat_L = latitude(latitude);
    sigma_hl_L = local_hour(lhour);
    sigma_day_L = day_year(dayear);
    sigma_alti_L = altitude(200.);
    sigma_FP_L = solar_flux(fpsolar);
    sigma_KP_L = magnetic_index(kpindex);
    stdev = sigma_alat_L + sigma_hl_L + sigma_day_L + sigma_alti_L +
            sigma_FP_L + sigma_KP_L + Copt_L - Res_L;
  } else if (altitude() <= 500e0) {
    /* Global 1-sigma calculation by means of statistical binning data (GOCE,
     * CHAMP, GRACE and SWARM)) and Hyperplan equation (6 dimensions)
     * Polynomial interpolation */
    Copt_M = 0.929e0;
    Res_M = 112.79e0;
    sigma_alat_M = latitude(latitude);
    sigma_hl_M = localHour(lhour);
    sigma_day_M = dayYear(dayear);
    sigma_alti_M = altitude(altitude);
    sigma_FP_M = solarFlux(fpsolar);
    sigma_KP_M = magneticIndex(kpindex);
    stdev = sigma_alat_M + sigma_hl_M + sigma_day_M + sigma_alti_M +
            sigma_FP_M + sigma_KP_M + Copt_M - Res_M;
  } else if (altitude() > 500e0) {
    /* Global 1-sigma calculation by means of statistical binning data (GOCE,
     * CHAMP, GRACE and SWARM) and Hyperplan equation (6 dimensions)
     * Polynomial extrapolation by fixing altitude value at 500 km */
    Copt_H = 0.929;
    Res_H = 112.79;
    sigma_alat_H = latitude(latitude);
    sigma_hl_H = localHour(lhour);
    sigma_day_H = dayYear(dayear);
    sigma_alti_H = altitude(500.);
    sigma_FP_H = solarFlux(fpsolar);
    sigma_KP_H = magneticIndex(kpindex);
    stdev = sigma_alat_H + sigma_hl_H + sigma_day_H + sigma_alti_H +
            sigma_FP_H + sigma_KP_H + Copt_H - Res_H;
  }
}
#endif
