#include "nrlmsise00.hpp"

using namespace dso::nrlmsise00::detail;

// if anything other than zero is returned, error!
int dso::Nrlmsise00::ghp7(const InParamsCore *in,
                          dso::nrlmsise00::OutParams *out,
                          double press) noexcept {

  constexpr const double bm = 1.3806e-19;

  double pl = std::log10(press);

  double zi = 0e0, z;
  // initial altitude estimate
  if (pl >= -5e0) {
    if (pl > 2.5e0)
      zi = 18.06e0 * (3.00e0 - pl);
    else if (pl > 0.75e0 && pl <= 2.5e0)
      zi = 14.98e0 * (3.08e0 - pl);
    else if (pl > -1.0e0 && pl <= 0.75e0)
      zi = 17.8e0 * (2.72e0 - pl);
    else if (pl > -2.0e0 && pl <= -1.0e0)
      zi = 14.28e0 * (3.64e0 - pl);
    else if (pl > -4.0e0 && pl <= -2.0e0)
      zi = 12.72e0 * (4.32e0 - pl);
    else if (pl <= -4.0e0)
      zi = 25.3 * (0.11e0 - pl);

    const double cl = in->glat / 90e0;
    const double cl2 = cl * cl;

    const double cd = (in->doy >= 182) ? ((in->doy / 91.25e0) - 3.0e0)
                                       : ((1e0 - in->doy) / 91.25e0);

    double ca = 0e0;
    if (pl > -0.23e0)
      ca = (2.79e0 - pl) / (2.79e0 + 0.23e0);
    else if (pl > -1.11e0 && pl <= -0.23e0)
      ca = 1e0;
    else if (pl <= -1.11e0 && pl > -3.0e0)
      ca = (-2.93e0 - pl) / (-2.93e0 + 1.11e0);

    z = zi - 4.87e0 * cl * cd * ca - 1.64e0 * cl2 * ca + 0.31e0 * ca * cl;
  }

  if (pl < -5e0)
    z = 22e0 * std::pow(pl + 4e0, 2e0) + 110e0;

  constexpr const int ltest = 12;
  constexpr const double test = 0.00043e0;
  const bool imr = in->meters();
  int l = 0;
  do {
    ++l;
    gtd7(in, out, 48);
    const double xn = out->d[0] + out->d[1] + out->d[2] + out->d[3] +
                      out->d[4] + out->d[6] + out->d[7];
    double p = bm * xn * out->t[1];
    if (imr)
      p *= 1e-6; //[m]
    const double diff = pl - std::log10(p);

    if (std::abs(diff) < test) {
      return 0;
    }

    if (l == ltest) {
      // Non converging, should not happen
      // alt = z;
      return 1;
    }

    double xm = out->d[5] / xn / 1.66e-24;
    if (in->meters())
      xm *= 1e3;
    const double g = gsurf / std::pow(1e0 + z / re, 2);
    const double sh = r100gas * out->t[1] / (xm * g);
    // New altitude estimate using scale height
    if (l < 6)
      z -= sh * diff * 2.302e0;
    else
      z -= sh * diff;
  } while (true);

  return 0;
}
