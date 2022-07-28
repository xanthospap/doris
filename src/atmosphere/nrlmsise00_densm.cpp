#include "nrlmsise00.hpp"
#include <algorithm>
#ifdef dEBUG
#include <cassert>
#endif

using namespace dso::nrlmsise00;

/// @brief Calculate Temperature and Density Profiles for lower atmos.
double dso::Nrlmsise00::densm(double alt, double d0, double xm,
                              double &tz) const noexcept {
  constexpr const int dim = 10;
  double xs[dim], ys[dim], y2out[dim], work[dim];

  double densm = d0;
  if (alt <= zn2[0]) {
    // STRATOSPHERE/MESOSPHERE TEMPERATURE
    double z = std::max(alt, zn2[mn2 - 1]);
    int mn = mn2;
    double z1 = zn2[0];
    double z2 = zn2[mn - 1];
    double t1 = tn2[0];
    double t2 = tn2[mn - 1];
    double zg = zeta(z, z1);
    double zgdif = zeta(z2, z1);

// setup spline nodes
#ifdef DEBUG
    assert(mn < 10);
#endif
    for (int k = 0; k < mn; k++) {
      xs[k] = zeta(zn2[k], z1) / zgdif;
    }
    for (int k = 0; k < mn; k++) {
      ys[k] = 1e0 / tn2[k];
    }
    double yd1 = -tgn2[0] / (t1 * t1) * zgdif;
    double yd2 =
        -tgn2[1] / (t2 * t2) * zgdif * std::pow((re + z2) / (re + z1), 2e0);

    // calculate spline coefficients
    dso::nrlmsise00::spline(xs, ys, mn, yd1, yd2, y2out, work);
    double x = zg / zgdif;
    double y = dso::nrlmsise00::splint(xs, ys, y2out, mn, x);

    // temperature at altitude
    tz = 1e0 / y;

    if (std::abs(xm) > nearzero) {
      // CALCULATE STRATOSPHERE/MESOSPHERE DENSITY
      const double glb = gsurf / std::pow(1e0 + z1 / re, 2e0);
      const double gamm = xm * glb * zgdif / r100gas;

      // integrate temperature profile
      const double yi = dso::nrlmsise00::splini(xs, ys, y2out, mn, x);
      const double expl = (gamm * yi > 50e0) ? (50e0) : (gamm * yi);

      // density at altitude
      densm = densm * (t1 / tz) * std::exp(-expl);
    } // if (std::abs(xm) > nearzero)

    if (alt <= zn3[0]) {
      // ! TROPOSPHERE/STRATOSPHERE TEMPERATURE
      z = alt;
      mn = mn3;
      z1 = zn3[0];
      z2 = zn3[mn - 1];
      t1 = tn3[0];
      t2 = tn3[mn - 1];
      zg = zeta(z, z1);
      zgdif = zeta(z2, z1);

      // Setup spline nodes
      for (int k = 0; k < mn; k++)
        xs[k] = zeta(zn3[k], z1) / zgdif;
      for (int k = 0; k < mn; k++)
        ys[k] = 1e0 / tn3[k];
      yd1 = -tgn3[0] / (t1 * t1) * zgdif;
      yd2 = -tgn3[1] / (t2 * t2) * zgdif * std::pow((re + z2) / (re + z1), 2e0);

      // calculate spline coefficients
      dso::nrlmsise00::spline(xs, ys, mn, yd1, yd2, y2out, work);
      x = zg / zgdif;
      y = dso::nrlmsise00::splint(xs, ys, y2out, mn, x);

      // temperature at altitude
      tz = 1e0 / y;
      if (std::abs(xm) > nearzero) {
        // CALCULATE TROPOSPHERIC/STRATOSPHERE DENSITY
        const double glb = gsurf / std::pow(1e0 + z1 / re, 2e0);
        const double gamm = xm * glb * zgdif / r100gas;

        // integrate temperature profile
        const double yi = dso::nrlmsise00::splini(xs, ys, y2out, mn, x);
        const double expl = (gamm * yi > 50e0) ? (50e0) : (gamm * yi);

        // density at altitude
        densm = densm * (t1 / tz) * std::exp(-expl);
      } // if (std::abs(xm) > nearzero)
    }   // if (alt <= zn3[0])
  }     // if (alt <= zn2[0])

  if (std::abs(xm) < nearzero)
    densm = tz;

  return densm;
}
