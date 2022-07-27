#include "geodesy/units.hpp"
#include "nrlmsise00.hpp"
#include <cstring>
#include <cstdio> // only for debuging

using namespace dso::nrlmsise00;

const double altl[8] = {200e0, 300e0, 160e0, 250e0,
                          240e0, 450e0, 320e0, 450e0};
const double alpha[9] = {-0.380e0, 0e0, 0e0, 
0e0, 0.170e0, 0e0, 
-0.380e0, 0e0, 0e0};
const int mt[11] = {48, 0, 4, 16, 28, 32, 40, 1, 49, 14, 17};
double zn1[5] = {120e0, 110e0, 100e0, 90e0, 72.5e0};


int dso::Nrlmsise00::gts7(const InParams *in, int mass,
                          OutParams *out) noexcept {
  
  // last call with this altitude
  static double altlast = -999e0;

  // input always changed from previous call ....
  constexpr const bool input_changed = true;
  
  printf("gts7 called ...\n");

  const double za = pdl[1][15];
  zn1[0] = za;
  std::memset(out->d, 0, sizeof(double) * 9);

  // TINF variations not important below za or zn1(1)
  double tinf = 0e0;
  if (in->alt > zn1[0]) {
    if (input_changed || altlast <= zn1[0]) {
      //printf("* calling globe7 from (1)");
      tinf = ptm[0] * pt[0] * (1e0 + in->sw.sw[15] * globe7(in, pt));
    }
  } else {
    tinf = ptm[0] * pt[0];
  }
  out->t[0] = tinf;
  //printf("t(1) = %+.15e\n", out->t[0]);

  // gradient variations not important below zn1(5)
  if (in->alt > zn1[4]) {
    if (input_changed || altlast <= zn1[4]) {
      //printf("* calling globe7 from (2)");
      xg0 = ptm[3] * ps[0] * (1e0 + in->sw.sw[18] * globe7(in, ps));
    }
  } else {
    xg0 = ptm[3] * ps[0];
  }

  // calculate these temperatures only if input changed
  //printf("-> checks : %1d and %10.2f\n", input_changed, in->alt);
  if (input_changed || in->alt < 300e0) {
    //printf("* calling globe7 from (3)");
    tlb = ptm[1] * (1e0 + in->sw.sw[16] * globe7(in, pd[3])) * pd[3][0];
    //printf("Computed tlb=%+25.14e\n", tlb);
  }
  const double s = xg0 / (tinf - tlb);

  // Lower thermosphere temp variations not significant for
  // density above 300 km
  if (in->alt >= 300e0) {
    tn1[1] = ptm[6] * ptl[0][0];
    tn1[2] = ptm[2] * ptl[1][0];
    tn1[3] = ptm[7] * ptl[2][0];
    tn1[4] = ptm[4] * ptl[3][0];
    tgn1[1] = ptm[8] * pma[8][0] * tn1[4] * tn1[4] /
              std::pow(ptm[4] * ptl[3][0], 2e0);
  } else if (input_changed || altlast >= 300e0) {
    tn1[1] = ptm[6] * ptl[0][0] / (1e0 - in->sw.sw[17] * glob7s(in, ptl[0]));
    tn1[2] = ptm[2] * ptl[1][0] / (1e0 - in->sw.sw[17] * glob7s(in, ptl[1]));
    tn1[3] = ptm[7] * ptl[2][0] / (1e0 - in->sw.sw[17] * glob7s(in, ptl[2]));
    tn1[4] = ptm[4] * ptl[3][0] /
             (1e0 - in->sw.sw[17] * in->sw.sw[19] * glob7s(in, ptl[3]));
    tgn1[1] = ptm[8] * pma[8][0] *
              (1e0 + in->sw.sw[17] * in->sw.sw[19] * glob7s(in, pma[8])) *
              tn1[4] * tn1[4] / std::pow(ptm[4] * ptl[3][0], 2);
  }

  [[maybe_unused]] const double z0 = zn1[3];
  [[maybe_unused]] const double t0 = tn1[3];
  [[maybe_unused]] const double tr12 = 1e0;

  int j = 0;
  if (mass) {

    // N2 variation factor at Zlb
    const double g28 = in->sw.sw[20] * globe7(in, pd[2]);

    // variation of turbopause height
    const double zhf =
        pdl[1][24] *
        (1e0 + in->sw.sw[4] * pdl[0][24] * std::sin(dso::deg2rad(in->glat)) *
                   std::cos(dr * (in->doy - pt[13])));
    out->t[0] = tinf;
    const double xmm = pdm[2][4];
    const double z = in->alt;

    bool goto_100 = false;
    int i;
    for (i = 0; i < 11; i++) {
      if (mt[i] == mass) {
        j = i+1;
        goto_100 = true;
        printf("j=%i\n", j);
        break;
      }
    }

    double tz = 0, zhm28 = 0, b28 = 0;
    double db01, db04, db14, db16, db28, db32, db40;
    if (goto_100) {
      printf("Computing N2 density d(3)\n");
      if (z < altl[5] || mass == 28 || mass == 48) {
        //
        //  **** N2 DENSITY ****
        //
        // diffusive density at zlb
        // printf("db28: %25.15f %25.15f %25.15f\n", pdm[2][0], g28, pd[2][0]);
        db28 = pdm[2][0] * std::exp(g28) * pd[2][0];
        // diffusive density at alt
        //printf("densu args: %15.8f %25.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n", 
        //z, db28, tinf, tlb, 28e0, alpha[2], out->t[1], ptm[5], s);
        out->d[2] =
            densu(z, db28, tinf, tlb, 28e0, alpha[2], out->t[1], ptm[5], s, zn1);
        printf("N2 density=%25.15f", out->d[2]);
        dd = out->d[2];
        // turbopause
        const double zh28 = pdm[2][2] * zhf;
        zhm28 = pdm[2][3] * pdl[1][5];
        const double xmd = 28e0 - xmm;
        // mixed density at Zlb
        //printf("densu args: %15.8f %25.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n", 
        //zh28, db28, tinf, tlb, xmd, alpha[2] - 1e0, tz, ptm[5], s);
        b28 = densu(zh28, db28, tinf, tlb, xmd, alpha[2] - 1e0, tz, ptm[5], s, zn1);
        //printf("b28=%25.15f\n", b28);
        if (z <= altl[2] && std::abs(in->sw.sw[14]) > 0e0) {
          // mixed density at alt
          dm28 = densu(z, b28, tinf, tlb, xmm, alpha[2], tz, ptm[5], s, zn1);
          // net density at alt
          out->d[2] = dnet(out->d[2], dm28, zhm28, xmm, 28e0);
          printf(" new %25.15f", out->d[2]);
        }
      }
      printf("\n");

      if (j==1 || j == 3 || j == 4 || j == 9) {
        printf("Computing HE density d(1)\n");
        //
        // **** HE DENSITY ****
        // BP: --
        // Density variation factor at Zlb
        const double g4 = in->sw.sw[20] * globe7(in, pd[0]);
        // diffusive density at zlb
        db04 = pdm[0][0] * std::exp(g4) * pd[0][0];
        //printf("densu args: %15.8f %25.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n", 
        //z, db04, tinf, tlb, 4e0, alpha[0], out->t[1], ptm[5], s);
        out->d[0] =
            densu(z, db04, tinf, tlb, 4e0, alpha[0], out->t[1], ptm[5], s, zn1);
        printf("HE density=%25.15f", out->d[0]);
        dd = out->d[0];
        if (z <= altl[0] && std::abs(in->sw.sw[14]) > 0e0) {
          // turbopause
          const double zh04 = pdm[0][2];
          // mixed density at Zlb
          const double b04 = densu(zh04, db04, tinf, tlb, 4e0 - xmm,
                                   alpha[0] - 1e0, out->t[1], ptm[5], s, zn1);
          // mixed density at alt
          const double dm04 =
              densu(z, b04, tinf, tlb, xmm, 0e0, out->t[1], ptm[5], s, zn1);
          const double zhm04 = zhm28;
          // net density at alt
          out->d[0] = dnet(out->d[0], dm04, zhm04, xmm, 4e0);
          // correction to specified mixing ratio at ground
          const double rl = std::log(b28 * pdm[0][1] / b04);
          const double zc04 = pdm[0][4] * pdl[1][0];
          const double hc04 = pdm[0][5] * pdl[1][1];
          // net density corrected at alt
          out->d[0] *= ccor(z, rl, hc04, zc04);
          printf(" new %25.15f", out->d[0]);
        }
        printf("\n");

        //
        // **** O DENSITY ****
        //
        // BP: --
        // Density variation factor at Zlb
        printf("Computing O density d(2)\n");
        const double g16 = in->sw.sw[20] * globe7(in, pd[1]);
        // diffusion density at Zlb
        db16 = pdm[1][0] * std::exp(g16) * pd[1][0];
        // diffusive density at alt
        out->d[1] =
            densu(z, db16, tinf, tlb, 16e0, alpha[1], out->t[1], ptm[5], s, zn1);
        printf(" O density=%25.15f", out->d[1]);
        dd = out->d[1];
        if (z <= altl[1] && std::abs(in->sw.sw[14]) > 0) {
          // corrected pdm(31) to pdm(3,2) 12/2/85
          // turbopause
          const double zh16 = pdm[1][2];
          // mixed density at Zlb
          const double b16 = densu(zh16, db16, tinf, tlb, 16e0 - xmm,
                                   alpha[1] - 1e0, out->t[1], ptm[5], s, zn1);
          // mixed density at alt
          const double dm16 =
              densu(z, b16, tinf, tlb, xmm, 0e0, out->t[1], ptm[5], s, zn1);
          const double zhm16 = zhm28;
          // net density at alt
          out->d[1] = dnet(out->d[1], dm16, zhm16, xmm, 16e0);
          // 3/16/99 Change form to match O2 departure from diff equil near
          // 150 km and add dependence on F10.7
          // RL=ALOG(B28*PDM(2,2)*ABS(PDL(17,2))/B16)
          const double rl =
              pdm[1][1] * pdl[1][16] *
              (1e0 + in->sw.sw[0] * pdl[0][23] * (in->f107A - 150e0));
          const double hc16 = pdm[1][5] * pdl[1][3];
          const double zc16 = pdm[1][4] * pdl[1][2];
          const double hc216 = pdm[1][5] * pdl[1][4];
          out->d[1] *= ccor2(z, rl, hc16, zc16, hc216);
          // chemistry correction
          const double hcc16 = pdm[1][7] * pdl[1][13];
          const double zcc16 = pdm[1][6] * pdl[1][12];
          const double rc16 = pdm[1][3] * pdl[1][14];
          // net density corrected at alt
          out->d[1] *= ccor(z, rc16, hcc16, zcc16);
          printf(" new %25.15f", out->d[1]);
        }
        printf("\n");
      }

      if (j==1 || j == 3 || j == 4 || j == 6 || j == 9) {

        //
        //  **** O2 DENSITY ****
        //
        printf("Computing O2 density\n");
        // BP: 200
        const double g32 = in->sw.sw[20] * globe7(in, pd[4]);
        // diffusive density at Zlb
        db32 = pdm[3][0] * std::exp(g32) * pd[4][0];
        // diffusive density at alt
        out->d[3] =
            densu(z, db32, tinf, tlb, 32e0, alpha[3], out->t[1], ptm[5], s, zn1);
        if (mass == 49) {
          dd += 2e0 * out->d[3];
        } else {
          dd = out->d[3];
        }
        if (std::abs(in->sw.sw[14]) > 0e0) {
          if (z <= altl[3]) {
            // turbopause
            const double zh32 = pdm[3][2];
            // mixed density at Zlb
            const double b32 = densu(zh32, db32, tinf, tlb, 32e0 - xmm,
                                     alpha[3] - 1e0, out->t[1], ptm[5], s, zn1);
            // mixed density at alt
            const double dm32 =
                densu(z, b32, tinf, tlb, xmm, 0e0, out->t[1], ptm[5], s, zn1);
            const double zhm32 = zhm28;
            // net density at alt
            out->d[3] = dnet(out->d[3], dm32, zhm32, xmm, 32e0);
            // correction to specified mixing ration at ground
            const double rl = std::log(b28 * pdm[3][1] / b32);
            const double hc32 = pdm[3][5] * pdl[1][7];
            const double zc32 = pdm[3][4] * pdl[1][6];
            out->d[3] *= ccor(z, rl, hc32, zc32);
          }
          // correction for general departure from diffusive
          // equilibrium above Zlb
          const double hcc32 = pdm[3][7] * pdl[1][22];
          const double hcc232 = pdm[3][7] * pdl[0][22];
          const double zcc32 = pdm[3][6] * pdl[1][21];
          const double rc32 =
              pdm[3][3] * pdl[1][23] *
              (1e0 + in->sw.sw[0] * pdl[0][23] * (in->f107A - 150e0));
          // net density corrected at alt
          out->d[3] *= ccor2(z, rc32, hcc32, zcc32, hcc232);
        }
      }
    }

    if (j==1 || j == 3 || j == 4 || j == 6 || j == 9 || j == 7) {

      //
      // **** AR DENSITY ****
      //
      printf("Computing AR density\n");
      // BP: 300
      const double g40 = in->sw.sw[20] * globe7(in, pd[5]);
      // diffusive density at Zlb
      db40 = pdm[4][0] * std::exp(g40) * pd[5][0];
      // diffusive density at alt
      out->d[4] =
          densu(z, db40, tinf, tlb, 40e0, alpha[4], out->t[1], ptm[5], s, zn1);
      dd = out->d[4];
      if (z <= altl[4] && std::abs(in->sw.sw[14]) > 0e0) {
        // turbopause
        const double zh40 = pdm[4][2];
        // mixed density at Zlb
        const double b40 = densu(zh40, db40, tinf, tlb, 40e0 - xmm,
                                 alpha[4] - 1e0, out->t[1], ptm[5], s, zn1);
        // mixed density at alt
        const double dm40 =
            densu(z, b40, tinf, tlb, xmm, 0e0, out->t[1], ptm[5], s, zn1);
        const double zhm40 = zhm28;
        // net density at alt
        out->d[4] = dnet(out->d[4], dm40, zhm40, xmm, 40e0);
        // correction to specified mixing ratio at ground
        const double rl = std::log(b28 * pdm[4][1] / b40);
        const double hc40 = pdm[4][5] * pdl[1][9];
        const double zc40 = pdm[4][4] * pdl[1][8];
        // net density corrected at alt
        out->d[4] *= ccor(z, rl, hc40, zc40);
      }
    }

    if (j==1 || j == 3 || j == 4 || j == 6 || j == 9 || j == 7 || j == 8) {
      //
      // **** HYDROGEN DENSITY ****
      printf("Computing Hy density\n");
      //
      // BP: 400
      const double g1 = in->sw.sw[20] * globe7(in, pd[6]);
      // diffusive density at zlb
      db01 = pdm[5][0] * std::exp(g1) * pd[6][0];
      // diffusive density at alt
      out->d[6] =
          densu(z, db01, tinf, tlb, 1e0, alpha[6], out->t[1], ptm[5], s, zn1);
      dd = out->d[6];
      if (z <= altl[6] && std::abs(in->sw.sw[14]) > 0e0) {
        // turbopause
        const double zh01 = pdm[5][2];
        // mixed density at Zlb
        const double b01 = densu(zh01, db01, tinf, tlb, 1e0 - xmm,
                                 alpha[6] - 1e0, out->t[1], ptm[5], s, zn1);
        // mixed density at alt
        const double dm01 =
            densu(z, b01, tinf, tlb, xmm, 0e0, out->t[1], ptm[5], s, zn1);
        const double zhm01 = zhm28;
        // net density at alt
        out->d[6] = dnet(out->d[6], dm01, zhm01, xmm, 1e0);
        // correction to specified mixing ratio at ground
        const double rl =
            std::log(b28 * pdm[5][1] * std::abs(pdl[1][17]) / b01);
        const double hc01 = pdm[5][5] * pdl[1][11];
        const double zc01 = pdm[5][4] * pdl[1][10];
        out->d[6] *= ccor(z, rl, hc01, zc01);
        // chemistry correction
        const double hcc01 = pdm[5][7] * pdl[1][19];
        const double zcc01 = pdm[5][6] * pdl[1][18];
        const double rc01 = pdm[5][3] * pdl[1][20];
        // net density corrected at alt
        out->d[6] *= ccor(z, rc01, hcc01, zcc01);
      }
    }

    if (j==1 || j == 3 || j == 4 || j == 6 || j == 9 || j == 7 || j == 8) {

      //
      // **** ATOMIC NITROGEN DENSITY ****
      printf("Computing Ni density d(8)\n");
      //
      // BP: 500
      const double g14 = in->sw.sw[20] * globe7(in, pd[7]);
      // diffusive density at Zlb
      db14 = pdm[6][0] * std::exp(g14) * pd[7][0];
      // diffusive density at alt
      out->d[7] =
          densu(z, db14, tinf, tlb, 14e0, alpha[7], out->t[1], ptm[5], s, zn1);
      printf("Ni: %25.12f", out->d[7]);
      dd = out->d[7];
      if (z <= altl[7] && std::abs(in->sw.sw[14]) > 0e0) {
        // turbopause
        const double zh14 = pdm[6][2];
        // mixed density at Zlb
        const double b14 = densu(zh14, db14, tinf, tlb, 14e0 - xmm,
                                 alpha[7] - 1e0, out->t[1], ptm[5], s, zn1);
        // mixed density at alt
        const double dm14 =
            densu(z, b14, tinf, tlb, xmm, 0e0, out->t[1], ptm[5], s, zn1);
        const double zhm14 = zhm28;
        // net density at alt
        out->d[7] = dnet(out->d[7], dm14, zhm14, xmm, 14e0);
        printf(" new %25.12f", out->d[7]);
        // correction to specified mixing ratio at ground
        //printf("rl: %25.15f %25.15f %25.15f %25.15f\n",b28, pdm[6][1], pdl[0][2],b14);
        const double rl = std::log(b28 * pdm[6][1] * std::abs(pdl[0][2]) / b14);
        const double hc14 = pdm[6][5] * pdl[0][1];
        const double zc14 = pdm[6][4] * pdl[0][0];
        //printf("Calling ccor: %25.15f %25.15f %25.15f %25.15f\n", z, rl, hc14, zc14);
        out->d[7] *= ccor(z, rl, hc14, zc14);
        printf(" new %25.12f", out->d[7]);
        // chemistry correction
        const double hcc14 = pdm[6][7] * pdl[0][4];
        const double zcc14 = pdm[6][6] * pdl[0][3];
        const double rc14 = pdm[6][3] * pdl[0][5];
        // net density corrected at alt
        out->d[7] *= ccor(z, rc14, hcc14, zcc14);
        printf(" new %25.12f", out->d[7]);
      }
      printf("\n");
    }
    if (j==1 || j == 3 || j == 4 || j == 6 || j == 9 || j == 7 || j == 8 || j == 10 ||
        j == 11) {

      //
      // **** Anomalous OXYGEN DENSITY ****
      //
      // BP: 600
      printf("Computing AO density\n");
      const double g16h = in->sw.sw[20] * globe7(in, pd[8]);
      const double db16h = pdm[7][0] * std::exp(g16h) * pd[8][0];
      const double tho = pdm[7][9] * pdl[0][6];
      dd = densu(z, db16h, tho, tho, 16e0, alpha[8], out->t[1], ptm[5], s, zn1);
      const double zsht = pdm[7][5];
      const double zmho = pdm[7][4];
      const double zsho = scalh(zmho, 16e0, tho);
      out->d[8] =
          dd * std::exp(-zsht / zsho * (std::exp(-(z - zmho) / zsht) - 1e0));
      if (mass == 48) {
        // total mass density
        out->d[5] =
            1.66e-24 * (4e0 * out->d[0] + 16e0 * out->d[1] + 28e0 * out->d[2] +
                        32e0 * out->d[3] + 40e0 * out->d[4] + out->d[6] +
                        14e0 * out->d[7]);
        [[maybe_unused]] const double db48 =
            1.66e-24 * (4e0 * db04 + 16e0 * db16 + 28e0 * db28 + 32e0 * db32 +
                        40e0 * db40 + db01 + 14e0 * db14);
      }
    }
  } // if mass

  [[maybe_unused]] double ddum;
  if (j != 5) {
    // BP: 700
    const double z = std::abs(in->alt);
    ddum = densu(z, 1e0, tinf, tlb, 0e0, 0e0, out->t[1], ptm[5], s, zn1);
  }

  if (in->meters()) {
    for (int i = 0; i < 9; i++) {
      out->d[i] *= 1e6;
    }
    out->d[5] /= 1e3;
  }
  altlast = in->alt;
  
  //for (int i=0; i<9; i++)
  //printf("d(%i) = %25.11e\n", i+1, out->d[i]);

  return 0;
}
