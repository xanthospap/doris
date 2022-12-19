#include "nrlmsise00.hpp"
#include <algorithm>
#include <cstring>

using namespace dso::nrlmsise00::detail;

int dso::Nrlmsise00::gtd7(const InParamsCore *in,
                          dso::nrlmsise00::OutParams *out, int mass) noexcept {
  static double altlast = 99999e0;
  constexpr const double zmix = 62.5e0;

  dso::nrlmsise00::OutParams outc;
  double *__restrict__ d = out->d;
  double *__restrict__ t = out->t;
  double *__restrict__ ds = outc.d;

  // test for changed input; fuck that, assert input changed, why give the same
  // input ?
  constexpr const bool input_changed = true;

  // latitude variation of gravity (none for SW(2)=0)
  double xlat = in->glat;
  if (std::abs(in->sw.sw[1]) < 0)
    xlat = 45e0;
  re = glatf(xlat, gsurf);

  const double xmm = pdm[2][4];

  // THERMOSPHERE/MESOSPHERE (above ZN2(1))
  const double altt = std::max(in->alt, zn2[0]);
  int mss = mass;
  // Only calculate N2 in thermosphere if alt in mixed region
  if (in->alt < zmix && mass > 0)
    mss = 28;
  // Only calculate thermosphere if input parameters changed (but they always
  // have in this C++ version ...) or altitude above ZN2(1) in mesosphere
  {
    InParamsCore inc(*in);
    inc.alt = altt;
    gts7(&inc, &outc, mss);
    dm28m = dm28;
    if (in->meters())
      dm28m = dm28 * 1e6;
    t[0] = outc.t[0];
    t[1] = outc.t[1];
    if (in->alt >= zn2[0]) {
      std::memcpy(d, outc.d, sizeof(double) * 9);
      return 0;
    }
  }

  // LOWER MESOSPHERE/UPPER STRATOSPHERE [between ZN3(1) and ZN2(1)]
  // Temperature at nodes and gradients at end nodes
  // Inverse temperature a linear function of spherical harmonics
  // Only calculate nodes if input changed
  const double sw19 = in->sw.sw[19];
  const double sw21 = in->sw.sw[21];
  if (input_changed || altlast >= zn2[0]) {
    tgn2[0] = tgn1[1];
    tn2[0] = tn1[4];
    tn2[1] = pma[0][0] * pavgm[0] / (1e0 - sw19 * glob7s(in, pma[0]));
    tn2[2] = pma[1][0] * pavgm[1] / (1e0 - sw19 * glob7s(in, pma[1]));
    tn2[3] = pma[2][0] * pavgm[2] / (1e0 - sw19 * sw21 * glob7s(in, pma[2]));
    tgn2[1] = pavgm[8] * pma[9][0] * (1e0 + sw19 * sw21 * glob7s(in, pma[9])) *
              tn2[3] * tn2[3] / std::pow(pma[2][0] * pavgm[2], 2e0);
  }

  // LOWER STRATOSPHERE AND TROPOSPHERE [below ZN3(1)]
  // Temperature at nodes and gradients at end nodes
  // Inverse temperature a linear function of spherical harmonics
  // Only calculate nodes if input changed
  if (in->alt < zn3[0]) {
    if (input_changed || altlast >= zn3[0]) {
      tgn3[0] = tgn2[1];
      tn3[1] = pma[3][0] * pavgm[3] / (1e0 - sw21 * glob7s(in, pma[3]));
      tn3[2] = pma[4][0] * pavgm[4] / (1e0 - sw21 * glob7s(in, pma[4]));
      tn3[3] = pma[5][0] * pavgm[5] / (1e0 - sw21 * glob7s(in, pma[5]));
      tn3[4] = pma[6][0] * pavgm[6] / (1e0 - sw21 * glob7s(in, pma[6]));
      tgn3[1] = pma[7][0] * pavgm[7] * (1e0 + sw21 * glob7s(in, pma[7])) *
                tn3[4] * tn3[4] / std::pow(pma[6][0] * pavgm[6], 2e0);
    }
  }

  double dmr, tz;
  if (mass == 0) {
    dd = densm(in->alt, 1e0, 0e0, tz);
    t[1] = tz;
  } else {
    // LINEAR TRANSITION TO FULL MIXING BELOW ZN2(1)
    double dmc = 0e0;
    if (in->alt > zmix)
      dmc = 1e0 - (zn2[0] - in->alt) / (zn2[0] - zmix);
    double dz28 = ds[2];
    // N2 DENSITY
    dmr = ds[2] / dm28m - 1e0;
    d[2] = densm(in->alt, dm28m, xmm, tz);
    d[2] *= (1e0 + dmr * dmc);
    // HE DENSITY
    d[1] = 0e0;
    if (mass == 4 || mass == 48) {
      dmr = ds[0] / (dz28 * pdm[0][1]) - 1e0;
      d[0] = d[2] * pdm[0][1] * (1e0 + dmr * dmc);
    }
    // O DENSITY
    d[1] = 0e0;
    d[8] = 0e0;
    // O2 DENSITY
    d[3] = 0e0;
    if (mass == 32 || mass == 48) {
      dmr = ds[3] / (dz28 * pdm[3][1]) - 1e0;
      d[3] = d[2] * pdm[3][1] * (1e0 + dmr * dmc);
    }
    // AR DENSITY
    d[4] = 0e0;
    if (mass == 40 || mass == 48) {
      dmr = ds[4] / (dz28 * pdm[4][1]) - 1e0;
      d[4] = d[2] * pdm[4][1] * (1e0 + dmr * dmc);
    }
    // HYDROGEN DENSITY
    d[6] = 0e0;
    // ATOMIC NITROGEN DENSITY
    d[7] = 0e0;

    // TOTAL MASS DENSITY
    if (mass == 48) {
      d[5] = 1.66e-24 * (4e0 * d[0] + 16e0 * d[1] + 28e0 * d[2] + 32e0 * d[3] +
                         40e0 * d[4] + d[6] + 14e0 * d[7]);
      if (in->meters())
        d[5] /= 1e3;
    }
    t[1] = tz;
  }

  // last used altitude
  altlast = in->alt;

  return 0;
}
