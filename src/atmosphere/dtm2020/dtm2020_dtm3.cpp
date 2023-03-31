#include "atmosphere/dtm2020/dtm2020.hpp"
#include <cstring>
#include "geodesy/units.hpp"

namespace {
constexpr const double cpmg = .19081;
constexpr const double re = 6356.77;
constexpr const double rgas = 831.4;
constexpr const double zlb0 = 120.;
constexpr const double xlmg = -1.2392;
/*constexpr const double cose = .9175;*/
constexpr const double gsurf = 980.665;
/*constexpr const double sine = .3978;*/
constexpr const double spmg = .98163;
/*constexpr const double zero = 0.;*/
const double alefa[] = {-0.40, -0.38, 0., 0., 0., 0.};
const int ma[] = {1, 4, 16, 28, 32, 14};
const double vma[] = {1.6606e-24,  6.6423e-24,  26.569e-24,
                      46.4958e-24, 53.1381e-24, 23.2479e-24};
} // unnamed namespace

/*
 * This is a translated version of the FORTAN source code found at
 * https://github.com/swami-h2020-eu/mcm/blob/main/src/dtm2020/dtm2020_F107_Kp-subr.f90
 *
 * Original (FORTRAN) Author and translated Version:
 * Author : SB (CNES)
 * last version updated : 21/10/2020
 *
 * Roriginal documentation block follows:
 * ***********************************************************************
 *  CNES DTM2020 operational: F10.7 and Kp
 *  ver 21/10/2020
 *  calculation of temperature and density with DTM2020_Oper
 *
 * *par ** INPUT **
 *      day = day of year [1-366]
 *      f   = f(1)=instantaneous flux at (t - 24hr)    /   f(2)=0.
 *      fbar= fbar(1)=mean flux of last 81 days at t   /   fbar(2)=0.
 *      akp = akp(1)= kp delayed by 3 hours, akp(3)=mean of last 24 hours
 *            akp(2) & akp(4)=0.
 *      alti= altitude (in km) greater than 120 km
 *      hl  = local time (in radian: 0-24hr = 0-2pi)
 *      alat= latitude (in radian)
 *      xlon= longitude (in radian)
 *
 * *par ** OUTPUT **
 *      tz      = temperature at altitude -> alti
 *      tinf    = exospheric temperature
 *      d(1)    = partial density of atomic hydrogen (in gram/cm3)
 *      d(2)    = partial density of helium
 *      d(3)    = partial density of atomic oxygen
 *      d(4)    = partial density of molecular nitrogen
 *      d(5)    = partial density of molecular oxygen
 *      d(6)    = partial density of atomic nitrogen
 *      ro      = total density (in gram/cm3)
 *      wmm     = mean molecular mass (in gram)
 * ***********************************************************************
 */
int dso::Dtm2020::dtm3() noexcept {
  /* data pool */
  double RawMem[nlatm * 9];
  /* initialize to zero */
  std::memset(RawMem, 0, nlatm * 9 * sizeof(double));
  /* only access RawMem from array pointers, all arrays already seto to 0 */
  double *__restrict__ daz  = RawMem + 0 * nlatm;
  double *__restrict__ daz2 = RawMem + 1 * nlatm;
  double *__restrict__ dh   = RawMem + 2 * nlatm;
  double *__restrict__ dhe  = RawMem + 3 * nlatm;
  double *__restrict__ do2  = RawMem + 4 * nlatm;
  double *__restrict__ dox  = RawMem + 5 * nlatm;
  double *__restrict__ dt0  = RawMem + 6 * nlatm;
  double *__restrict__ dtp  = RawMem + 7 * nlatm;
  double *__restrict__ dtt  = RawMem + 8 * nlatm;

  const double zlb = zlb0;

  const double c = std::sin(in.latitude());
  const double c2 = c * c;
  const double c4 = c2 * c2;
  const double s = cos(in.latitude());
  const double s2 = s * s;

  /* following are instance variables */
  p10 = c;
  p20 = 1.5 * c2 - 0.5;
  p30 = c * (2.5 * c2 - 1.5);
  p40 = 4.375 * c4 - 3.75 * c2 + 0.375;
  p50 = c * (7.875 * c4 - 8.75 * c2 + 1.875);
  p60 = (5.5 * c * p50 - 2.5 * p40) / 3.;
  p11 = s;
  p21 = 3. * c * s;
  p31 = s * (7.5 * c2 - 1.5);
  p41 = c * s * (17.5 * c2 - 7.5);
  p51 = s * (39.375 * c4 - 26.25 * c2 + 1.875);
  p22 = 3. * s2;
  p32 = 15. * c * s2;
  p42 = s2 * (52.5 * c2 - 7.5);
  p52 = 3. * c * p42 - 2. * p32;
  p62 = 2.75 * c * p52 - 1.75 * p42;
  p33 = 15. * s * s2;

  const double clmlmg = std::cos(in.longitude() - xlmg);
  const double sp = s * cpmg * clmlmg + c * spmg;
  const double cmg = sp;
  const double cmg2 = cmg * cmg;
  const double cmg4 = cmg2 * cmg2;

  /* following are instance variables */
  p10mg = cmg;
  p20mg = 1.5 * cmg2 - 0.5;
  p40mg = 4.375 * cmg4 - 3.75 * cmg2 + 0.375;
  hl0 = in.local_hours_radians();
  ch = std::cos(hl0);
  sh = std::sin(hl0);
  c2h = ch*ch - sh*sh;
  s2h = 2. * ch * sh;
  c3h = c2h*ch - s2h*sh;
  s3h = s2h*ch + c2h*sh;

  const double gdelt = gldtm(tt, dtt, 1., in.longitude());
  dtt[0] = 1. + gdelt;
  out.exospheric_temperature() = tt[0] * dtt[0]; /* exospheric temperature */
  const double gdelt0 = gldtm(t0, dt0, 1., in.longitude());
  dt0[0] = 1. + gdelt0;
  const double t120 = t0[0] * dt0[0];
  const double gdeltp = gldtm(tp, dtp, 1., in.longitude());
  dtp[0] = 1. + gdeltp;
  const double tp120 = tp[0] * dtp[0];

  const double sigma = tp120 / (out.exospheric_temperature() - t120);
  const double dzeta = (re + zlb) / (re + in.altitude_km());
  const double zeta = (in.altitude_km() - zlb) * dzeta;
  /* const double dzeta2 = dzeta * dzeta; never used */
  const double sigzeta = sigma * zeta;
  const double expsz = std::exp(-sigzeta);
  /* temperature at altitude */
  out.temperature() = out.exospheric_temperature() - (out.exospheric_temperature() - t120) * expsz;

  double dbase[6];
  const double gdelh = gldtm(h, dh, 0., in.longitude());
  dh[0] = std::exp(gdelh);
  dbase[0] = h[0] * dh[0];

  const double gdelhe = gldtm(he, dhe, 0., in.longitude());
  dhe[1] = std::exp(gdelhe);
  dbase[1] = he[0] * dhe[0];

  const double gdelo = gldtm(o, dox, 1., in.longitude());
  dox[0] = std::exp(gdelo);
  dbase[2] = o[0] * dox[0];

  const double gdelaz2 = gldtm(az2, daz2, 1., in.longitude());
  daz2[0] = std::exp(gdelaz2);
  dbase[3] = az2[0] * daz2[0];

  const double gdelo2 = gldtm(o2, do2, 1., in.longitude());
  do2[0] = std::exp(gdelo2);
  dbase[4] = o2[0] * do2[0];

  const double gdelaz = gldtm(az, daz, 1., in.longitude());
  daz[0] = std::exp(gdelaz);
  dbase[5] = az[0] * daz[0];

  double glb = gsurf / std::pow(1. + zlb / re, 2e0);
  glb /= (sigma * rgas * out.exospheric_temperature());
  const double t120tz = t120 / out.temperature();
  /* following are never used */
  /* const double xlog = std::log(t120tz); */
  /* const double tinftz = out.exospheric_temperature() / out.temperature(); */
  /* const double t120tt = t120 / (out.exospheric_temperature() - t120);*/

  out.total_density() = 0e0;
  double cc[6];
  for (int i = 0; i < 6; i++) {
    const double gamma = ma[i] * glb;
    const double upapg = 1. + alefa[i] + gamma;
    const double fz = std::pow(t120tz,upapg) * std::exp(-sigzeta * gamma);
    cc[i] = dbase[i] * fz;
    const double d = cc[i] * vma[i]; /* partial density of ... */
    out.total_density() += d; /* total density */
  }
  
  /* mean molescular mass */
  out.mean_molecular_mass() =
      out.total_density() /
      (vma[0] * (cc[0] + cc[1] + cc[2] + cc[3] + cc[4] + cc[5]));

  return 0;
}
