#include "nrlmsise00.hpp"
#include <algorithm>
#include <cassert>
#include <cstdio> // only for debugging

using namespace dso::nrlmsise00;

/// @brief  Calculate Temperature and Density Profiles for MSIS models
/// New lower thermo polynomial 10/30/89
double dso::Nrlmsise00::densu(double alt, double dlb, double t1, double t2,
                              double xm, double xalph, double &tz, double zlb,
                              double s2, const double *zn1) noexcept {
  constexpr const int dim = 5;
  double xs[dim], ys[dim], y2out[dim], work[dim];

  double densu = 1e0;
  //printf("called densu ...\n");

  // joining altitude of Bates and spline
  const double za = zn1[0];
  double z = std::max(alt, za);

  // geopotential altitude difference from ZLB
  const double zg2 = zeta(z, zlb);

  // Bates temperature
  const double tt = t1 - (t1 - t2) * std::exp(-s2 * zg2);
  //printf("t1=%20.15e t2=%20.15e s2=%20.15e zg2=%20.15e\n", t1,t2,s2,zg2);
  const double ta = tt;
  tz = tt;
  //printf("setting tz=%+25.15e\n", tz);
  densu = tz;

  int mn;
  double z1, tt1, zgdif, x;
  if (alt < za) {
    // CALCULATE TEMPERATURE BELOW ZA
    // Temperature gradient at ZA from Bates profile
    const double dta = (t1 - ta) * s2 * std::pow((re + zlb) / (re + za), 2e0);
    tgn1[0] = dta;
    tn1[0] = ta;
    z = std::max(alt, zn1[mn1 - 1]);
    mn = mn1;
    z1 = zn1[0];
    const double z2 = zn1[mn - 1];
    tt1 = tn1[0];
    const double tt2 = tn1[mn - 1];

    // geopotential difference from Z1
    const double zg = zeta(z, z1);
    zgdif = zeta(z2, z1);
    //printf("zeta(%25.15f %25.15f)=%25.15f\n", z2,z1,zgdif);

    // setup spline nodes
    for (int k = 0; k < mn; k++)
      xs[k] = zeta(zn1[k], z1) / zgdif;
    for (int k = 0; k < mn; k++)
      ys[k] = 1e0 / tn1[k];

    // end node derivatives
    const double yd1 = -tgn1[0] / (tt1 * tt1) * zgdif;
    const double yd2 =
        -tgn1[1] / (tt2 * tt2) * zgdif * std::pow((re + z2) / (re + z1), 2e0);

    // calculate spline coefficients
#ifdef DEBUG
    assert(mn <= 5);
#endif
    dso::nrlmsise00::spline(xs, ys, mn, yd1, yd2, y2out, work);
    x = zg / zgdif;
    //printf("x=%25.15f / %25.15f\n", zg, zgdif);
    const double y = dso::nrlmsise00::splint(xs, ys, y2out, mn, x);

    // temperature at altitude
    tz = 1e0 / y;
    //printf("setting tz=%+25.15e\n", tz);
    densu = tz;
  } // if (alt < za)

  if (std::abs(xm) > nearzero) {
    // CALCULATE DENSITY ABOVE ZA
    double glb = gsurf / std::pow(1e0 + zlb / re, 2e0);
    const double gammo = xm * glb / (s2 * r100gas * t1);
    double expl = std::exp(-s2 * gammo * zg2);
    if (expl > 50e0 || tt <= 0e0)
      expl = 50e0;

    // density at altitude
    const double densa = dlb * std::pow(t2 / tt, 1e0 + xalph + gammo) * expl;
    //printf("setting densa=%+25.15e\n", densa);
    densu = densa;

    if (alt < za) {
      // CALCULATE DENSITY BELOW ZA
      glb = gsurf / std::pow(1e0 + z1 / re, 2e0);
      const double gamm = xm * glb * zgdif / r100gas;

      // integrate spline temperatures
      double yi = dso::nrlmsise00::splini(xs, ys, y2out, mn, x);
      expl = gamm * yi;
      if (expl > 50e0 || tz <= 0e0)
        expl = 50e0;

      // density at altitude
      densu = densu * std::pow(tt1 / tz, 1e0 + xalph) * std::exp(-expl);
      //printf("setting densu=%+25.15e\n", densu);
    } // if (alt<za)
  }   // if (std::abs(xm) > nearzero)

  return densu;
}