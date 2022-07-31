#include "geodesy/units.hpp"
#include "nrlmsise00.hpp"

using namespace dso::nrlmsise00::detail;

double dso::Nrlmsise00::glob7s(const InParamsCore *in, double *pp) noexcept {
  
  double t[14] = {0e0};
  const double glong = in->glon;
  const double doy = in->doy;

  // only access pp through here, let compiler know
  double *__restrict__ p = pp;

  // VERSION OF GLOBE FOR LOWER ATMOSPHERE 10/26/99
  constexpr const double pset = 2e0;
  static int last_doy = -1e0; // last doy used
  static double p32 = -1e3;
  static double p18 = -1e3;
  static double p14 = -1e3;
  static double p39 = -1e3;
  static double cd14;
  static double cd18;
  static double cd32;
  static double cd39;

  if (std::abs(p[99]) < nearzero)
    p[99] = pset;

  // did day of year change?
  const bool doy_changed = (int)in->doy - last_doy != 0;
  if (doy_changed || std::abs(p32 - p[31]) > nearzero)
    cd32 = std::cos(dr * (doy - p[31]));
  if (doy_changed || std::abs(p18 - p[17]) > nearzero)
    cd18 = std::cos(2e0 * dr * (doy - p[17]));
  if (doy_changed || std::abs(p14 - p[13]) > nearzero)
    cd14 = std::cos(dr * (doy - p[13]));
  if (doy_changed || std::abs(p39 - p[38]) > nearzero)
    cd39 = std::cos(2e0 * dr * (doy - p[38]));

  // update last used doy
  last_doy = doy;

  p32 = p[31];
  p18 = p[17];
  p14 = p[13];
  p39 = p[39];

  t[0] = p[21] * dfa;
  
  // time independent
  t[1] = p[1] * plg[0][2] + p[2] * plg[0][4] + p[22] * plg[0][6] +
         p[26] * plg[0][1] + p[14] * plg[0][3] + p[59] * plg[0][5];
  
  // symmetrical annual
  t[2] = (p[18] + p[47] * plg[0][2] + p[29] * plg[0][4]) * cd32;
  
  // symmetrical semi-annual
  t[3] = (p[15] + p[16] * plg[0][2] + p[30] * plg[0][4]) * cd18;
  
  // asymmetrical annual
  t[4] = (p[9] * plg[0][1] + p[10] * plg[0][3] + p[20] * plg[0][5]) * cd14;
  
  // asymmetric semi-annual
  t[5] = (p[37] * plg[0][1]) * cd39;

  double absw[Switches::dim];
  for (int i=0;i<Switches::dim; i++) absw[i] = std::abs(in->sw.sw[i]);

  // diurnal
  if (absw[6] > 0) {
    const double t71 = p[11] * plg[1][2] * cd14 * in->sw.swc[4];
    const double t72 = p[12] * plg[1][2] * cd14 * in->sw.swc[4];
    t[6] = ((p[3] * plg[1][1] + p[4] * plg[1][3] + t71) * ctloc +
            (p[6] * plg[1][1] + p[7] * plg[1][3] + t72) * stloc);
  }
  
  // semidiurnal
  if (absw[7] > 0) {
    const double t81 =
        (p[23] * plg[2][3] + p[35] * plg[2][5]) * cd14 * in->sw.swc[4];
    const double t82 =
        (p[33] * plg[2][3] + p[36] * plg[2][5]) * cd14 * in->sw.swc[4];
    t[7] = ((p[5] * plg[2][2] + p[41] * plg[2][4] + t81) * c2tloc +
            (p[8] * plg[2][2] + p[42] * plg[2][4] + t82) * s2tloc);
  }
  
  // terdiurnal
  if (absw[13] > 0) {
    t[13] = p[39] * plg[3][3] * s3tloc + p[40] * plg[3][3] * c3tloc;
  }
  
  // magnetic activity
  if (absw[8] > 0) {
    if (in->sw.sw[8] > 0)
      t[8] = apdf * (p[32] + p[45] * plg[0][2] * in->sw.swc[1]);
    if (in->sw.sw[8] < 0)
      t[8] = (p[50] * apt[0] + p[96] * plg[0][2] * apt[0] * in->sw.swc[1]);
  }
  
  // longitudinal
  if (absw[9] > 0 && absw[10] > 0 &&
      glong > -1e3) {
    t[10] = (1e0 +
             plg[0][1] *
                 (p[80] * in->sw.swc[4] * std::cos(dr * (doy - p[81])) +
                  p[85] * in->sw.swc[5] * std::cos(2e0 * dr * (doy - p[86]))) +
             p[83] * in->sw.swc[2] * std::cos(dr * (doy - p[84])) +
             p[87] * in->sw.swc[3] * std::cos(2e0 * dr * (doy - p[88]))) *
            ((p[64] * plg[1][2] + p[65] * plg[1][4] + p[66] * plg[1][6] +
              p[74] * plg[1][1] + p[75] * plg[1][3] + p[76] * plg[1][5]) *
                 std::cos(dso::deg2rad(glong)) +
             (p[90] * plg[1][2] + p[91] * plg[1][4] + p[92] * plg[1][6] +
              p[77] * plg[1][1] + p[78] * plg[1][3] + p[79] * plg[1][5]) *
                 std::sin(dso::deg2rad(glong)));
  }

  double tt = 0e0;
  for (int i = 0; i < 14; i++) {
    tt += absw[i] * t[i];
  }

  return tt;
}
