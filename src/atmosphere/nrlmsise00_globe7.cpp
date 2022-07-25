#include "geodesy/units.hpp"
#include "nrlmsise00.hpp"
#include "geodesy/units.hpp"

using namespace dso::nrlmsise00;

// TODO
// should be using doy here or fractional doy ? See line ~12

double dso::Nrlmsise00::globe7(const InParams *in, double *p) noexcept {
    
    // calculate G(L) function

    // upper thermosphere parameters
    constexpr const int nsw = 14;
    constexpr const double hr = 0.2618e0;
    constexpr const double sr = 7.2722e-5;

    static double xl = 1000e0;
    static double tll = 1000e0;
    static double sw9 = 1e0;
    static double dayl = -1e0;
    static double p14 = -1e3;
    static double p18 = -1e3;
    static double p32 = -1e3;
    static double p39 = -1e3;
    static double cd14, cd18, cd32, cd39;
    static double stloc, ctloc, s2tloc, c2tloc, s3tloc, c3tloc;

    double t[14] = {0e0};
    if (in->sw.sw[8] > 0e0) sw9 = 1e0;
    if (in->sw.sw[8] < 0e0) sw9 = -1e0;

    if (std::abs(xl-in->glat)>nearzero) {
        // calculate legendre polynomials
        const double c = std::sin(dso::deg2rad(in->glat));
        const double s = std::cos(dso::deg2rad(in->glat));
        const double c2 = c*c;
        const double c4 = c2*c2;
        const double s2 = s*s;
        plg[0][1] = c;
        plg[0][2] = 0.5e0  * (3e0 * c2 - 1e0);
        plg[0][3] = 0.5e0  * (5e0 * c * c2 - 3e0 * c);
        plg[0][4] = (35e0    * c4 - 30e0    * c2 + 3e0) / 8e0;
        plg[0][5] =
            (63e0    * c2 * c2 * c - 70e0    * c2 * c + 15e0    * c) / 8e0;
        plg[0][6] = (11e0    * c * plg[0][5] - 5e0 * plg[0][4]) / 6e0;
        // PLG(8,1) = (13.0_r8*C*PLG(7,1) - d_6*PLG(6,1))/7.0_r8
        plg[1][1] = s;
        plg[1][2] = 3e0 * c * s;
        plg[1][3] = 1.5e0  * (5e0 * c2 - 1e0) * s;
        plg[1][4] = 2.5e0  * (7e0 * c2 * c - 3e0 * c) * s;
        plg[1][5] = 1.875e0  * (21e0    * c4 - 14e0    * c2 + 1e0) * s;
        plg[1][6] = (11e0    * c * plg[1][5] - 6e0 * plg[1][4]) / 5e0;
        // PLG(8,2) = (13.0_r8*C*PLG(7,2)-d_7*PLG(6,2))/d_6
        // PLG(9,2) = (15.0_r8*C*PLG(8,2)-d_8*PLG(7,2))/d_7
        plg[2][2] = 3e0 * s2;
        plg[2][3] = 15e0    * s2 * c;
        plg[2][4] = 7.5e0  * (7e0 * c2 - 1e0) * s2;
        plg[2][5] = 3e0 * c * plg[2][4] - 2e0 * plg[2][3];
        plg[2][6] = (11e0    * c * plg[2][5] - 7e0 * plg[2][4]) / 4e0;
        plg[2][7] = (13e0    * c * plg[2][6] - 8e0 * plg[2][5]) / 5e0;
        plg[3][3] = 15e0    * s2 * s;
        plg[3][4] = 105e0    * s2 * s * c;
        plg[3][5] = (9e0 * c * plg[3][4] - 7e0 * plg[3][3]) / 2e0;
        plg[3][6] = (11e0    * c * plg[3][5] - 8e0 * plg[3][4]) / 3e0;
        xl = in->glat;
    }

    const double tloc = in->lst;
    if (std::abs(tll-tloc) > nearzero) {
        if (std::abs(in->sw.sw[6])>0e0 || std::abs(in->sw.sw[7])>0e0 || std::abs(in->sw.sw[13])) {
            stloc = std::sin(hr*tloc);
            ctloc = std::cos(hr*tloc);
            s2tloc = std::sin(2e0*hr*tloc);
            c2tloc = std::cos(2e0*hr*tloc);
            s3tloc = std::sin(3e0*hr*tloc);
            c3tloc = std::cos(3e0*hr*tloc);
            tll = tloc;
        }
    }

    const double day = in->doy;
    if (std::abs(day - dayl) > nearzero || std::abs(p[13] - p14) > nearzero)
      cd14 = std::cos(dr * (day - p[13]));
    if (std::abs(day - dayl) > nearzero || std::abs(p[17] - p18) > nearzero)
      cd18 = std::cos(2e0 * dr * (day - p[17]));
    if (std::abs(day - dayl) > nearzero || std::abs(p[31] - p32) > nearzero)
      cd32 = std::cos(dr * (day - p[31]));
    if (std::abs(day - dayl) > nearzero || std::abs(p[38] - p39) > nearzero)
      cd39 = std::cos(2e0 * dr * (day - p[38]));
    dayl = day;

    p14 = p[13];
    p18 = p[17];
    p32 = p[31];
    p39 = p[38];

    // F10.7 effect
    const double df = in->f107 - in->f107A;
    const double dfa = in->f107A - 150e0;
    t[0] = p[19] * df * (1e0 + p[59] * dfa) + p[20] * df * df + p[21] * dfa +
           p[29] * dfa * dfa;
    const double f1 =
        1e0 + (p[47] * dfa + p[19] * df + p[20] * df * df) * in->sw.swc[0];
    const double f2 =
        1e0 + (p[49] * dfa + p[19] * df + p[20] * df * df) * in->sw.swc[0];
    // time independent
    t[1] = (p[1] * plg[0][2] + p[2] * plg[0][4] + p[22] * plg[0][6]) +
           (p[14] * plg[0][2]) * dfa * in->sw.swc[0] + p[26] * plg[0][1];
    // symmetric annual
    t[2] = p[18] * cd32;
    // symmetric semiannual
    t[3] = (p[15] + p[16] * plg[0][2]) * cd18;
    // asymmetric annual
    t[4] = f1 * (p[9] * plg[0][1] + p[10] * plg[0][3]) * cd14;
    // asymmetric semiannual
    t[5] = p[37] * plg[0][1] * cd39;
    // diurnal
    if (std::abs(in->sw.sw[6]) > 0e0) {
      const double t71 = (p[11] * plg[1][2]) * cd14 * in->sw.swc[4];
      const double t72 = (p[12] * plg[1][2]) * cd14 * in->sw.swc[4];
      t[6] = f2 *
             ((p[3] * plg[1][1] + p[4] * plg[1][3] + p[27] * plg[1][5] + t71) *
                  ctloc +
              (p[6] * plg[1][1] + p[7] * plg[1][3] + p[28] * plg[1][5] + t72) *
                  stloc);
    }
    // semidiurnal
    if (std::abs(in->sw.sw[7]) > 0e0) {
      const double t81 = (p[23] * plg[2][3] + p[35] * plg[2][5]) * cd14 * in->sw.swc[4];
      const double t82 = (p[33] * plg[2][3] + p[36] * plg[2][5]) * cd14 * in->sw.swc[4];
      t[7] = f2 * ((p[5] * plg[2][2] + p[41] * plg[2][4] + t81) * c2tloc +
                   (p[8] * plg[2][2] + p[42] * plg[2][4] + t82) * s2tloc);
    }
    // terdiurnal
    if (std::abs(in->sw.sw[13]) > 0e0) {
      t[13] = f2 * ((p[39] * plg[3][3] +
                     (p[93] * plg[3][4] + p[46] * plg[3][6]) * cd14 * in->sw.swc[4]) *
                        s3tloc +
                    (p[40] * plg[3][3] +
                     (p[94] * plg[3][4] + p[48] * plg[3][6]) * cd14 * in->sw.swc[4]) *
                        c3tloc);
    }
    // magnetic activity based on daily ap
    if (std::abs(sw9 + 1e0)> nearzero) {
        const double apd = ap[0] - 4e0;
        double p44 = p[43];
        const double p45 = p[44];
        if (p44 < 0e0) p44 = 1e-5;
        const double apdf =
            apd + (p45 - 1e0) * (apd + (std::exp(-p44 * apd) - 1e0) / p44);
        if (std::abs(in->sw.sw[8])> nearzero) {
          t[8] =
              apdf *
              (p[32] + p[45] * plg[0][2] + p[34] * plg[0][4] +
               (p[100] * plg[0][1] + p[101] * plg[0][3] + p[102] * plg[0][5]) *
                   cd14 * in->sw.swc[4] +
               (p[121] * plg[1][1] + p[122] * plg[1][3] + p[123] * plg[1][5]) *
                   in->sw.swc[6] * std::cos(hr * (tloc - p[124])));
        }
    } else if (std::abs(p[51])>nearzero) {
      double exp1 =
          std::exp(-10800e0 * std::abs(p[51]) /
                   (1e0 + p[138] * (45e0 - std::abs(in->glat))));
      if (exp1 > 0.99999e0)
        exp1 = 0.99999;
        if (p[24]<1e-4) p[24] = 1e-4;
        apt[0] = sg0(exp1);
        // APT(2)=SG2(EXP1)
        // APT(3)=SG0(EXP2)
        // APT(4)=SG2(EXP2)
        if (std::abs(in->sw.sw[8])>0) {
          t[8] =
              apt[0] *
              (p[50] + p[96] * plg[0][2] + p[54] * plg[0][4] +
               (p[125] * plg[0][1] + p[126] * plg[0][3] + p[127] * plg[0][5]) *
                   cd14 * in->sw.swc[4] +
               (p[128] * plg[1][1] + p[129] * plg[1][3] + p[130] * plg[1][5]) *
                   in->sw.swc[6] * std::cos(hr * (tloc - p[131])));
        }
    }