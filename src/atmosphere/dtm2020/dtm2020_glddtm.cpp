#include "atmosphere/dtm2020/dtm2020.hpp"

/*
 * This is a translated version of the FORTAN source code found at
 * https://github.com/swami-h2020-eu/mcm/blob/main/src/dtm2020/dtm2020_F107_Kp-subr.f90
 * and
 * https://github.com/swami-h2020-eu/mcm/blob/main/src/dtm2020/dtm2020_F30_Hp-subr.f90
 *
 * Basically, the two version of the gldtm/gldtm_Hp functions defined in these
 * two files are pretty similar. Hence, i merged them!
 *
 * In the function defines here, the input (integer) variable `hp` can be set
 * to 1 to emulate the gldtm_Hp routine. Else, it will do what gldtm does.
 *
 * Original (FORTRAN) Author and translated Version:
 * Author :  (CNES)
 * last version updated : 21/10/2020
 */

namespace {
constexpr const double rot = .017214206e0;
constexpr const double rot2 = .034428412e0;
constexpr const int ikp = 62;
constexpr const int ikpm = 67;
/* the following do not appear to be needed */
/* constexpr const double roth = .261799387e0; */
/* constexpr const double rots = 7.27220e-05; */
} // unnamed namespace

/*
 * Fortran's function signatures are:
 * gldtm(f,fbar,akp,day,a,da,gdel,ff0,xlon)
 * gldtm_Hp(f,fbar,akp,day,a,da,gdel,ff0,xlon)
 *
 * Blocks commented out, are actually for the gldtm_Hp version of the code.
 *
 * @param[out] a   array of size nlatm
 * @param[out] da  array of size nlatm
 * @param[in] fbar array of size 2
 * @param[in] akp  array of size 4
 * @param[in] f    array of size 4
 * @param[in] ff0  scalar
 * @param[in] longitide scalar [-π, π] in [rad]
 *
 * @obsolete Hp   if set to 1, then the gldtm_Hp version is called; else
 *                (i.e. Hp==0) the gldtm version is called
 * @return gdel
 */
double dso::Dtm2020::gldtm(const double *const f, const double *const fbar,
                           const double *const akp,
                           const double *__restrict__ a,
                           double *__restrict__ da, double ff0,
                           double longitude) noexcept {
  double fmfb[] = {0e0, 0e0};
  fmfb[0] = f[0] - fbar[0];
  fmfb[1] = f[1] - fbar[1];

  double fbm150[] = {0e0, 0e0};
  fbm150[0] = fbar[0] - 150.;
  fbm150[1] = fbar[1];

  da[2 - 1] = p20;
  da[3 - 1] = p40;
  da[74 - 1] = p10;
  da[77 - 1] = p30;
  da[78 - 1] = p50;
  da[79 - 1] = p60;

  da[4 - 1] = fmfb[0];
  da[6 - 1] = fbm150[0];
  da[5 - 1] = da[4 - 1] * da[4 - 1];
  da[69 - 1] = da[6 - 1] * da[6 - 1];
  da[82 - 1] = da[4 - 1] * p10;
  da[83 - 1] = da[4 - 1] * p20;
  da[84 - 1] = da[4 - 1] * p30;
  da[85 - 1] = da[6 - 1] * p20;
  da[86 - 1] = da[6 - 1] * p30;
  da[87 - 1] = da[6 - 1] * p40;

  double dkp, dkpm, c2fi = 0e0;
  ;
  // if (Hp == 1) { /* Hp version here! */
  //   if (akp[1 - 1] >= 9e0 && akp[3 - 1] > 7.5e0)
  //     akp[1 - 1] = 9e0 + (akp[1 - 1] - 9e0) / 5e0;
  //   if (akp[4 - 1] >= 9e0 && akp[3 - 1] > 7.5e0)
  //     akp[4 - 1] = 9e0 + (akp[4 - 1] - 9e0) / 5e0;
  //   if (akp[5 - 1] >= 9e0 && akp[3 - 1] > 7.5e0)
  //     akp[5 - 1] = 9e0 + (akp[5 - 1] - 9e0) / 5e0;
  //   if (akp[6 - 1] >= 9e0 && akp[3 - 1] > 7.5e0)
  //     akp[6 - 1] = 9e0 + (akp[6 - 1] - 9e0) / 5e0;
  //   if (akp[7 - 1] >= 9e0 && akp[3 - 1] > 7.5e0)
  //     akp[7 - 1] = 9e0 + (akp[7 - 1] - 9e0) / 5e0;
  //   if (akp[8 - 1] >= 9e0 && akp[3 - 1] > 7.5e0)
  //     akp[8 - 1] = 9e0 + (akp[8 - 1] - 9e0) / 5e0;
  //   dkp = akp[0];
  //   dkpm = akp[2];
  // } else { /* non-Hp version */
  c2fi = 1e0 - p10mg * p10mg;
  dkp = akp[1 - 1] + (a[ikp - 1] + c2fi * a[ikp]) * akp[2 - 1];
  const double dakp =
      a[7 - 1] + a[8 - 1] * p20mg + a[68 - 1] * p40mg +
      2e0 * dkp * (a[60 - 1] + a[61 - 1] * p20mg + a[75 - 1] * 2e0 * dkp * dkp);
  da[ikp - 1] = dakp * akp[2 - 1];
  da[ikp] = da[ikp - 1] * c2fi;
  dkpm = akp[3 - 1] + a[ikpm - 1] * akp[4 - 1];
  const double dakpm =
      a[64 - 1] + a[65 - 1] * p20mg + a[72 - 1] * p40mg +
      2e0 * dkpm *
          (a[66 - 1] + a[73 - 1] * p20mg + a[76 - 1] * 2e0 * dkpm * dkpm);
  da[ikpm - 1] = dakpm * akp[4 - 1];
  //}

  da[7 - 1] = dkp;
  da[8 - 1] = p20mg * dkp;
  da[60 - 1] = dkp * dkp;
  da[61 - 1] = p20mg * da[60 - 1];
  da[64 - 1] = dkpm;
  da[65 - 1] = p20mg * dkpm;
  da[66 - 1] = dkpm * dkpm;
  da[68 - 1] = p40mg * dkp;

  // if (Hp == 1) { /* Hp version */
  //   da[62 - 1] = p30mg * dkp;
  //   da[63 - 1] = p10mg * dkp;
  //   da[67 - 1] = p60mg * dkp;
  //   da[73 - 1] = p20mg * da[66 - 1];
  //   const double flux = std::max(f[0], fbar[0]);
  //   int iflux = static_cast<int>(flux);
  //   if (iflux >= 200) {
  //     da[75 - 1] = 0.333 * da[60 - 1] * da[60 - 1];
  //     da[76 - 1] = 0.1 * akp[7 - 1] * akp[7 - 1] * akp[7 - 1] * akp[7 - 1];
  //     da[79 - 1] = 0.1 * akp[8 - 1] * akp[8 - 1] * akp[8 - 1] * akp[8 - 1];
  //     da[71 - 1] = 0.1 * akp[5 - 1] * akp[5 - 1] * akp[5 - 1] * akp[5 - 1];
  //     da[72 - 1] = 0.1 * akp[6 - 1] * akp[6 - 1] * akp[6 - 1] * akp[6 - 1];
  //   } else if (iflux >= 190 && iflux < 200) {
  //     da[75 - 1] = 0.55 * da[60 - 1] * da[60 - 1];
  //     da[76 - 1] = 0.15 * akp[7 - 1] * akp[7 - 1] * akp[7 - 1] * akp[7 - 1];
  //     da[79 - 1] = 0.15 * akp[8 - 1] * akp[8 - 1] * akp[8 - 1] * akp[8 - 1];
  //     da[71 - 1] = 0.15 * akp[5 - 1] * akp[5 - 1] * akp[5 - 1] * akp[5 - 1];
  //     da[72 - 1] = 0.15 * akp[6 - 1] * akp[6 - 1] * akp[6 - 1] * akp[6 - 1];
  //   } else if (iflux >= 180 && iflux < 190) {
  //     da[75 - 1] = 0.733 * da[60 - 1] * da[60 - 1];
  //     da[76 - 1] = 0.2 * akp[7 - 1] * akp[7 - 1] * akp[7 - 1] * akp[7 - 1];
  //     da[79 - 1] = 0.2 * akp[8 - 1] * akp[8 - 1] * akp[8 - 1] * akp[8 - 1];
  //     da[71 - 1] = 0.2 * akp[5 - 1] * akp[5 - 1] * akp[5 - 1] * akp[5 - 1];
  //     da[72 - 1] = 0.2 * akp[6 - 1] * akp[6 - 1] * akp[6 - 1] * akp[6 - 1];
  //   } else if (iflux >= 160 && iflux < 180) {
  //     da[75 - 1] = da[60 - 1] * da[60 - 1];
  //     da[76 - 1] = 0.4 * akp[7 - 1] * akp[7 - 1] * akp[7 - 1] * akp[7 - 1];
  //     da[79 - 1] = 0.4 * akp[8 - 1] * akp[8 - 1] * akp[8 - 1] * akp[8 - 1];
  //     da[71 - 1] = 0.4 * akp[5 - 1] * akp[5 - 1] * akp[5 - 1] * akp[5 - 1];
  //     da[72 - 1] = 0.4 * akp[6 - 1] * akp[6 - 1] * akp[6 - 1] * akp[6 - 1];
  //   } else if (iflux >= 140 && iflux < 160) {
  //     da[75 - 1] = da[60 - 1] * da[60 - 1];
  //     da[76 - 1] = 0.8 * akp[7 - 1] * akp[7 - 1] * akp[7 - 1] * akp[7 - 1];
  //     da[79 - 1] = 0.8 * akp[8 - 1] * akp[8 - 1] * akp[8 - 1] * akp[8 - 1];
  //     da[71 - 1] = 0.8 * akp[5 - 1] * akp[5 - 1] * akp[5 - 1] * akp[5 - 1];
  //     da[72 - 1] = 0.8 * akp[6 - 1] * akp[6 - 1] * akp[6 - 1] * akp[6 - 1];
  //   } else {
  //     da[75 - 1] = da[60 - 1] * da[60 - 1];
  //     da[71 - 1] = akp[5 - 1] * akp[5 - 1] * akp[5 - 1] * akp[5 - 1];
  //     da[72 - 1] = akp[6 - 1] * akp[6 - 1] * akp[6 - 1] * akp[6 - 1];
  //     da[76 - 1] = akp[7 - 1] * akp[7 - 1] * akp[7 - 1] * akp[7 - 1];
  //     da[79 - 1] = akp[8 - 1] * akp[8 - 1] * akp[8 - 1] * akp[8 - 1];
  //   }
  //   da[70 - 1] = akp[2 - 1];
  // } else { /* non-Hp version */
  da[75 - 1] = da[60 - 1] * da[60 - 1];
  da[72 - 1] = p40mg * dkpm;
  da[73 - 1] = p20mg * da[66 - 1];
  da[76 - 1] = da[66 - 1] * da[66 - 1];
  //}

  double f0 = a[4 - 1] * da[4 - 1] + a[5 - 1] * da[5 - 1] +
              a[6 - 1] * da[6 - 1] + a[69 - 1] * da[69 - 1] +
              a[82 - 1] * da[82 - 1] + a[83 - 1] * da[83 - 1] +
              a[84 - 1] * da[84 - 1] + a[85 - 1] * da[85 - 1] +
              a[86 - 1] * da[86 - 1] + a[87 - 1] * da[87 - 1];
  const double f1f = 1e0 + f0 * ff0;
  f0 +=
      a[2 - 1] * da[2 - 1] + a[3 - 1] * da[3 - 1] + a[74 - 1] * da[74 - 1] +
      a[77 - 1] * da[77 - 1] + a[7 - 1] * da[7 - 1] + a[8 - 1] * da[8 - 1] +
      a[60 - 1] * da[60 - 1] + a[61 - 1] * da[61 - 1] + a[68 - 1] * da[68 - 1] +
      a[64 - 1] * da[64 - 1] + a[65 - 1] * da[65 - 1] + a[66 - 1] * da[66 - 1] +
      a[72 - 1] * da[72 - 1] + a[73 - 1] * da[73 - 1] + a[75 - 1] * da[75 - 1] +
      a[76 - 1] * da[76 - 1] + a[78 - 1] * da[78 - 1] + a[79 - 1] * da[79 - 1];
  /* in Hp case, add more ... */
  // if (Hp)
  //   f0 += (a[70 - 1] * da[70 - 1] + a[71 - 1] * da[71 - 1] +
  //          a[62 - 1] * da[62 - 1] + a[63 - 1] * da[63 - 1] +
  //          a[67 - 1] * da[67 - 1]);

  da[9 - 1] = std::cos(rot * (in.day_of_year() - a[11 - 1]));
  da[10 - 1] = p20 * da[9 - 1];
  da[12 - 1] = std::cos(rot2 * (in.day_of_year() - a[14 - 1]));
  da[13 - 1] = p20 * da[12 - 1];
  const double coste = std::cos(rot * (in.day_of_year() - a[18 - 1]));
  da[15 - 1] = p10 * coste;
  da[16 - 1] = p30 * coste;
  da[17 - 1] = da[6 - 1] * da[15 - 1];
  const double cos2te = std::cos(rot2 * (in.day_of_year() - a[20 - 1]));
  da[19 - 1] = p10 * cos2te;
  da[39 - 1] = p30 * cos2te;
  da[59 - 1] = da[6 - 1] * da[19 - 1];
  da[21 - 1] = p11 * ch;
  da[22 - 1] = p31 * ch;
  da[23 - 1] = da[6 - 1] * da[21 - 1];
  da[24 - 1] = da[21 - 1] * coste;
  da[25 - 1] = p21 * ch * coste;
  da[26 - 1] = p11 * sh;
  da[27 - 1] = p31 * sh;
  da[28 - 1] = da[6 - 1] * da[26 - 1];
  da[29 - 1] = da[26 - 1] * coste;
  da[30 - 1] = p21 * sh * coste;
  da[94 - 1] = p51 * ch;
  da[95 - 1] = p51 * sh;
  da[31 - 1] = p22 * c2h;
  da[37 - 1] = p42 * c2h;
  da[32 - 1] = p32 * c2h * coste;
  da[33 - 1] = p22 * s2h;
  da[38 - 1] = p42 * s2h;
  da[34 - 1] = p32 * s2h * coste;
  da[88 - 1] = p32 * c2h;
  da[89 - 1] = p32 * s2h;
  da[90 - 1] = da[6 - 1] * da[31 - 1];
  da[91 - 1] = da[6 - 1] * da[33 - 1];
  da[92 - 1] = p62 * c2h;
  da[93 - 1] = p62 * s2h;
  da[35 - 1] = p33 * c3h;
  da[36 - 1] = p33 * s3h;
  double fp =
      a[9 - 1] * da[9 - 1] + a[10 - 1] * da[10 - 1] + a[12 - 1] * da[12 - 1] +
      a[13 - 1] * da[13 - 1] + a[15 - 1] * da[15 - 1] + a[16 - 1] * da[16 - 1] +
      a[17 - 1] * da[17 - 1] + a[19 - 1] * da[19 - 1] + a[21 - 1] * da[21 - 1] +
      a[22 - 1] * da[22 - 1] + a[23 - 1] * da[23 - 1] + a[24 - 1] * da[24 - 1] +
      a[25 - 1] * da[25 - 1] + a[26 - 1] * da[26 - 1] + a[27 - 1] * da[27 - 1] +
      a[28 - 1] * da[28 - 1] + a[29 - 1] * da[29 - 1] + a[30 - 1] * da[30 - 1] +
      a[31 - 1] * da[31 - 1] + a[32 - 1] * da[32 - 1] + a[33 - 1] * da[33 - 1] +
      a[34 - 1] * da[34 - 1] + a[35 - 1] * da[35 - 1] + a[36 - 1] * da[36 - 1] +
      a[37 - 1] * da[37 - 1] + a[38 - 1] * da[38 - 1] + a[39 - 1] * da[39 - 1] +
      a[59 - 1] * da[59 - 1] + a[88 - 1] * da[88 - 1] + a[89 - 1] * da[89 - 1] +
      a[90 - 1] * da[90 - 1] + a[91 - 1] * da[91 - 1] + a[92 - 1] * da[92 - 1] +
      a[93 - 1] * da[93 - 1] + a[94 - 1] * da[94 - 1] + a[95 - 1] * da[95];
  da[40 - 1] = p10 * coste * dkp;
  da[41 - 1] = p30 * coste * dkp;
  da[42 - 1] = p50 * coste * dkp;
  da[43 - 1] = p11 * ch * dkp;
  da[44 - 1] = p31 * ch * dkp;
  da[45 - 1] = p51 * ch * dkp;
  da[46 - 1] = p11 * sh * dkp;
  da[47 - 1] = p31 * sh * dkp;
  da[48 - 1] = p51 * sh * dkp;
  fp +=
      a[40 - 1] * da[40 - 1] + a[41 - 1] * da[41 - 1] + a[42 - 1] * da[42 - 1] +
      a[43 - 1] * da[43 - 1] + a[44 - 1] * da[44 - 1] + a[45 - 1] * da[45 - 1] +
      a[46 - 1] * da[46 - 1] + a[47 - 1] * da[47 - 1] + a[48 - 1] * da[48 - 1];
  // if (!Hp) {
  //   const double dakp =
  //       (a[40 - 1] * p10 + a[41 - 1] * p30 + a[42 - 1] * p50) * coste +
  //       (a[43 - 1] * p11 + a[44 - 1] * p31 + a[45 - 1] * p51) * ch +
  //       (a[46 - 1] * p11 + a[47 - 1] * p31 + a[48 - 1] * p51) * sh;
  //   da[ikp - 1] = da[ikp - 1] + dakp * akp[1];
  //   da[ikp] = da[ikp - 1] + dakp * c2fi * akp[1];
  // }
  const double clfl = std::cos(longitude);
  da[49 - 1] = p11 * clfl;
  da[50 - 1] = p21 * clfl;
  da[51 - 1] = p31 * clfl;
  da[52 - 1] = p41 * clfl;
  da[53 - 1] = p51 * clfl;
  const double slfl = std::sin(longitude);
  da[54 - 1] = p11 * slfl;
  da[55 - 1] = p21 * slfl;
  da[56 - 1] = p31 * slfl;
  da[57 - 1] = p41 * slfl;
  da[58 - 1] = p51 * slfl;

  fp += a[49 - 1] * da[49 - 1] + a[50 - 1] * da[50 - 1] +
        a[51 - 1] * da[51 - 1] + a[52 - 1] * da[52 - 1] +
        a[53 - 1] * da[53 - 1] + a[54 - 1] * da[54 - 1] +
        a[55 - 1] * da[55 - 1] + a[56 - 1] * da[56 - 1] +
        a[57 - 1] * da[57 - 1] + a[58 - 1] * da[58];

  /* result */
  const double gdel = f0 + fp * f1f;
  return gdel;
}
