#include <cmath>
#include <algorithm>
#include "iers2010/iersc.hpp"

struct Jb2008InParams {
    double f10, f10b, s10, s10b, xm10, xm10b, y10, y10b, dstdtc;
};

struct Jb2008OutParams {
    double temp[2];
    double rho;
};

double xambar(double z) noexcept {

  const double c[] = {28.15204e0, -8.5586e-2, +1.2840e-4, -1.0056e-5,
                         -1.0210e-5, +1.5044e-6, +9.9826e-8};

  // Evaluates Equation (1)
  const double dz = z - 100e0;
  double amb = c[6];
  for (int i=6; i>=0; i--) {
    amb = dz * amb + c[i];
  }
  return amb;
}

double xgrav(double z) noexcept {
   return 9.80665e0/std::pow(1e0 + z/6356.766e0, 2e0);
}

double xlocal(double z, const double *tc) noexcept {
  const double dz = z - 125e0;
  double xl = 0e0;
  if (dz > 0e0) {
    xl=
        tc[0] + tc[2] * datan(tc[3] * dz * (1e0 + 4.5e-6 * std::pow(dz, 2.5)));
  } else {
    xl =
        ((-9.8204695e-6 * dz - 7.3039742e-4) * dz * dz + 1e0) * dz * tc[1] +
        tc[0];
  }

  return xl;
}

/// @brief Compute dTc correction for Jacchia-Bowman model
/// @param[in]  f10   f10 flux
/// @param[in]  xlst  local solar time (hours 0-23.999)
/// @param[in]  xlat  xlat = sat lat (rad)
/// @param[in]  zht   zht = height (km)
/// @param[out] dtc   dtc correction
void dtsub(double f10, double xlst, double xlat, double zht,
           double &dtc) noexcept {
  const double b[] = {
      -0.457512297e+01, -0.512114909e+01, -0.693003609e+02, 0.203716701e+03,
      0.703316291e+03,  -0.194349234e+04, 0.110651308e+04,  -0.174378996e+03,
      0.188594601e+04,  -0.709371517e+04, 0.922454523e+04,  -0.384508073e+04,
      -0.645841789e+01, 0.409703319e+02,  -0.482006560e+03, 0.181870931e+04,
      -0.237389204e+04, 0.996703815e+03,  0.361416936e+02};
  const double c[] = {
      -0.155986211e+02, -0.512114909e+01, -0.693003609e+02, 0.203716701e+03,
      0.703316291e+03,  -0.194349234e+04, 0.110651308e+04,  -0.220835117e+03,
      0.143256989e+04,  -0.318481844e+04, 0.328981513e+04,  -0.135332119e+04,
      0.199956489e+02,  -0.127093998e+02, 0.212825156e+02,  -0.275555432e+01,
      0.110234982e+02,  0.148881951e+03,  -0.751640284e+03, 0.637876542e+03,
      0.127093998e+02,  -0.212825156e+02, 0.275555432e+01};

  dtc = 0e0;

  const double tx = xlst / 24e0;
  const double ycs = std::cos(xlat);
  const double f = (f10 - 100e0) / 100e0;

  // calculates dTc
  const double tx2 = tx * tx;
  const double tx3 = tx2 * tx;
  const double tx4 = tx3 * tx;
  const double tx5 = tx4 * tx;

  if (zht >= 120e0 && zht <= 200e0) {
    const double h = (zht - 200e0) / 50e0;
    const double dtc200 = +c[16] + c[17] * tx * ycs + c[18] * tx2 * ycs +
                          c[19] * tx3 * ycs + c[20] * f * ycs +
                          c[21] * tx * f * ycs + c[22] * tx2 * f * ycs;
    const double sum =
        c[0] + b[1] * f + c[2] * tx * f + c[3] * tx2 * f + c[4] * tx3 * f +
        c[5] * tx4 * f + c[6] * tx5 * f + c[7] * tx * ycs + c[8] * tx2 * ycs +
        c[9] * tx3 * ycs + c[10] * tx4 * ycs + c[11] * tx5 * ycs + c[12] * ycs +
        c[13] * f * ycs + c[14] * tx * f * ycs + c[15] * tx2 * f * ycs;
    const double dtc200dz = sum;
    const double cc = 3e0 * dtc200 - dtc200dz;
    const double dd = dtc200 - cc;
    const double zp = (zht - 120e0) / 80e0;
    dtc = cc * zp * zp + dd * zp * zp * zp;
  
  } else if (zht > 200e0 && zht <= 240e0) {
    
    const double h = (zht - 200e0) / 50e0;
    const double sum =
        c[0] * h + b[1] * f * h + c[2] * tx * f * h + c[3] * tx2 * f * h +
        c[4] * tx3 * f * h + c[5] * tx4 * f * h + c[6] * tx5 * f * h +
        c[7] * tx * ycs * h + c[8] * tx2 * ycs * h + c[9] * tx3 * ycs * h +
        c[10] * tx4 * ycs * h + c[11] * tx5 * ycs * h + c[12] * ycs * h +
        c[13] * f * ycs * h + c[14] * tx * f * ycs * h +
        c[15] * tx2 * f * ycs * h + c[16] + c[17] * tx * ycs +
        c[18] * tx2 * ycs + c[19] * tx3 * ycs + c[20] * f * ycs +
        c[21] * tx * f * ycs + c[22] * tx2 * f * ycs;
    dtc = sum;
  
  } else if (zht > 240e0 && zht <= 300e0) {
    
    double h = (40e0) / 50e0;
    const double aa =
        c[0] * h + b[1] * f * h + c[2] * tx * f * h + c[3] * tx2 * f * h +
        c[4] * tx3 * f * h + c[5] * tx4 * f * h + c[6] * tx5 * f * h +
        c[7] * tx * ycs * h + c[8] * tx2 * ycs * h + c[9] * tx3 * ycs * h +
        c[10] * tx4 * ycs * h + c[11] * tx5 * ycs * h + c[12] * ycs * h +
        c[13] * f * ycs * h + c[14] * tx * f * ycs * h +
        c[15] * tx2 * f * ycs * h + c[16] + c[17] * tx * ycs +
        c[18] * tx2 * ycs + c[19] * tx3 * ycs + c[20] * f * ycs +
        c[21] * tx * f * ycs + c[22] * tx2 * f * ycs;
    const double bb =
        c[0] + b[1] * f + c[2] * tx * f + c[3] * tx2 * f + c[4] * tx3 * f +
        c[5] * tx4 * f + c[6] * tx5 * f + c[7] * tx * ycs + c[8] * tx2 * ycs +
        c[9] * tx3 * ycs + c[10] * tx4 * ycs + c[11] * tx5 * ycs + c[12] * ycs +
        c[13] * f * ycs + c[14] * tx * f * ycs + c[15] * tx2 * f * ycs;
    h = 300e0 / 100e0;
    const double dtc300 =
        b[0] + b[1] * f + b[2] * tx * f + b[3] * tx2 * f + b[4] * tx3 * f +
        b[5] * tx4 * f + b[6] * tx5 * f + b[7] * tx * ycs + b[8] * tx2 * ycs +
        b[9] * tx3 * ycs + b[10] * tx4 * ycs + b[11] * tx5 * ycs +
        b[12] * h * ycs + b[13] * tx * h * ycs + b[14] * tx2 * h * ycs +
        b[15] * tx3 * h * ycs + b[16] * tx4 * h * ycs + b[17] * tx5 * h * ycs +
        b[18] * ycs;
    const double dtc300dz = b[12] * ycs + b[13] * tx * ycs + b[14] * tx2 * ycs +
                            b[15] * tx3 * ycs + b[16] * tx4 * ycs +
                            b[17] * tx5 * ycs;
    const double cc = 3e0 * dtc300 - dtc300dz - 3e0 * aa - 2e0 * bb;
    const double dd = dtc300 - aa - bb - cc;
    const double zp = (zht - 240e0) / 60e0;
    dtc = aa + bb * zp + cc * zp * zp + dd * zp * zp * zp;
  
  } else if (zht > 300e0 && zht <= 600e0) {
    
    const double h = zht / 100e0;
    const double sum = b[0] + b[1] * f + b[2] * tx * f + b[3] * tx2 * f +
                       b[4] * tx3 * f + b[5] * tx4 * f + b[6] * tx5 * f +
                       b[7] * tx * ycs + b[8] * tx2 * ycs + b[9] * tx3 * ycs +
                       b[10] * tx4 * ycs + b[11] * tx5 * ycs + b[12] * h * ycs +
                       b[13] * tx * h * ycs + b[14] * tx2 * h * ycs +
                       b[15] * tx3 * h * ycs + b[16] * tx4 * h * ycs +
                       b[17] * tx5 * h * ycs + b[18] * ycs;
    dtc = sum;
  
  } else if (zht > 600e0 && zht <= 800e0) {
    
    const double zp = (zht - 600e0) / 100e0;
    const double hp = 600e0 / 100e0;
    const double aa = b[0] + b[1] * f + b[2] * tx * f + b[3] * tx2 * f +
                      b[4] * tx3 * f + b[5] * tx4 * f + b[6] * tx5 * f +
                      b[7] * tx * ycs + b[8] * tx2 * ycs + b[9] * tx3 * ycs +
                      b[10] * tx4 * ycs + b[11] * tx5 * ycs + b[12] * hp * ycs +
                      b[13] * tx * hp * ycs + b[14] * tx2 * hp * ycs +
                      b[15] * tx3 * hp * ycs + b[16] * tx4 * hp * ycs +
                      b[17] * tx5 * hp * ycs + b[18] * ycs;
    const double bb = b[12] * ycs + b[13] * tx * ycs + b[14] * tx2 * ycs +
                      b[15] * tx3 * ycs + b[16] * tx4 * ycs + b[17] * tx5 * ycs;
    const double cc = -(3e0 * aa + 4e0 * bb) / 4e0;
    const double dd = (aa + bb) / 4e0;
    dtc = aa + bb * zp + cc * zp * zp + dd * zp * zp * zp;
  }

  return;
}

/// Compute semiannual variation (delta log rho)
/// input day, height, f10b, s10b, m10b fsmb
///       025.  650.   150.  148.  147. 151.
/// output functions fz, gt, and del log rho value
///
/// @param[in]  day     day of year
/// @param[in]  ht      height (km)
/// @param[in]  f10b    ave 81-day centered f10
/// @param[in]  s10b    ave 81-day centered s10
/// @param[in]  xm10b   ave 81-day centered m10
/// @param[out] fzz     semiannual amplitude
/// @param[out] gtz     semiannual phase function
/// @param[out] drlog   delta log rho
void semian08(double day, double ht, double f10b, double s10b, double xm10b,
              double &fzz, double &gtz, double &drlog) noexcept {

  const double htz = ht / 1000e0;

  // compute new 81-day centered solar index for fz
  double fsmb = 1e0 * f10b - 0.70e0 * s10b - 0.04e0 * xm10b;

  // FZ blobal model values
  // 1997-2006 fit
  const double fzm[] = {0.2689e+00, -0.1176e-01, 0.2782e-01, -0.2782e-01,
                        0.3470e-03};

  fzz = fzm[0] + fzm[1] * fsmb + fzm[2] * fsmb * htz +
        fzm[3] * fsmb * htz * htz + fzm[4] * fsmb * fsmb * htz;

  if (fzz < 1e-6)
    fzz = 1e-6;

  // compute daily 81-day centered solar index for gt
  fsmb = 1e00 * f10b - 0.75e0 * s10b - 0.37e0 * xm10b;

  const double tau = (day - 1e0) / 365e0;
  const double sin1p = std::sin(TWOPI * tau);
  const double cos1p = std::cos(TWOPI * tau);
  const double sin2p = std::sin(2. * TWOPI * tau);
  const double cos2p = std::cos(2. * TWOPI * tau);

  // GT global model values
  // 1997-2006 fit:
  const double gtm[] = {-0.3633e+00, 0.8506e-01,  0.2401e+00, -0.1897e+00,
                        -0.2554e+00, -0.1790e-01, 0.5650e-03, -0.6407e-03,
                        -0.3418e-02, -0.1252e-02};

  gtz = gtm[0] + gtm[1] * sin1p + gtm[2] * cos1p + gtm[3] * sin2p +
        gtm[4] * cos2p + gtm[5] * fsmb + gtm[6] * fsmb * sin1p +
        gtm[7] * fsmb * cos1p + gtm[8] * fsmb * sin2p + gtm[9] * fsmb * cos2p;

  drlog = fzz * gtz;

  return;
}

/// @brief Compute day and year from time d1950 (days since 1950)
void tmoutd(double d1950, double &iyt, double &day) noexcept {
  double iyday = d1950;
  const double frac = d1950 - iyday;
  iyday = iyday + 364;
  double itemp = iyday / 1461;
  iyday = iyday - itemp * 1461;
  double iyr = 1949 + 4 * itemp;
  itemp = iyday / 365;
  if (itemp >= 3)
    itemp = 3;
  iyr = iyr + itemp;
  iyday = iyday - 365 * itemp + 1;
  iyr = iyr - 1900;
  day = iyday + frac;
  if (iyr >= 100)
    iyr = iyr - 100;
}

// The alpha are the thermal diffusion coefficients in Eq. (6)
const double alpha[] = {0e0,0e0,0e0,0e0,-0.38e0};

// AL10 is DLOG(10.0)
constexpr const double al10 = 2.3025851e0;

// The AMW are the molecular weights in order: N2, O2, O, Ar, He & H
const double amw[] = {28.0134e0, 31.9988e0, 15.9994e0,
                      39.9480e0, 4.0026e0,  1.00797e0};

// AVOGAD is Avogadro's number in mks units (molecules/kmol)
constexpr const double avogad = 02257e26;

// π related consts
constexpr const double PI = iers2010::DPI;
constexpr const double TWOPI = 2e0*iers2010::DPI;
constexpr const double FOURPI = 4e0*iers2010::DPI;
constexpr const double PIOV2 = 2e0/iers2010::DPI;
constexpr const double PIOV4 = 4e0/iers2010::DPI;
constexpr const double DEGRAD = PI / 180e0;

// The FRAC are the assumed sea-level volume fractions in order:
// N2, O2, Ar, and He
const double frac[] = {0.78110e0,0.20955e0,9.3400e-3,1.2890e-5};

// RSTAR is the universal gas-constant in mks units (joules/K/kmol)
constexpr const double rstar = 8314.32e0;

// The R# are values used to establish height step sizes in
// the regimes 90km to 105km, 105km to 500km and 500km upward.
constexpr const double r1 =0.010e0;
constexpr const double r2 =0.025e0;
constexpr const double r3 =0.075e0;

// The WT are weights for the Newton-Cotes Five-Point Quad. formula
const double wt[] = {0.311111111111111e0, 1.422222222222222e0,
                     0.533333333333333e0, 1.422222222222222e0,
                     0.311111111111111e0};

// The CHT are coefficients for high altitude density correction
const double cht[] = {0.22e0,-0.20e-02,0.115e-02,-0.211e-05};


int jb2008(double amjd, const double *sun, const double *sat,
           const Jb2008InParams &in, Jb2008OutParams &out) noexcept {

    double al10n[6], aln[6], tc[4];

    // Equation 14
    double fn = std::pow(in.f10b / 240.e0, 1e0 / 4e0);
    if (fn > 1e0)
      fn = 1e0;
    const double fsb = in.f10b * fn + in.s10b * (1e0 - fn);
    const double tsubc =
        392.4e0 + 3.227e0 * fsb + 0.298e0 * (in.f10 - in.f10b) +
        2.259e0 * (in.s10 - in.s10b) + 0.312e0 * (in.xm10 - in.xm10b) +
        0.178e0 * (in.y10 - in.y10b);
    
    // Equation 15
    const double eta =   0.5e0 * std::abs(sat[1] - sun[1]);
    const double theta = 0.5e0 * std::abs(sat[1] + sun[1]);

    // Equation 16
    const double h = sat[0] - sun[0];
    const double tau = h - 0.64577182e0 + 0.10471976e0 * std::sin(h + 0.75049158e0);
    const double glat  = sat[1];
    const double zht   = sat[2];
    const double glst  = h + PI;
    double glsthr = (glst/DEGRAD)*(24e0/360e0);
    if (glsthr >= 24e0) glsthr = glsthr - 24e0;
    if (glsthr < 0e0) glsthr = glsthr + 24e0;

    // Equation 17
    const double c = std::pow(std::cos(eta),2.5);
    const double s = std::pow(std::sin(theta),2.5);
    const double df = s + (c - s) * std::pow(std::abs(std::cos(0.5e0 * tau)),3e0);
    const double tsubl = tsubc * (1e0 + 0.31e0 * df);

    //  Compute correction to dTc for local solar time and lat correction
    dtsub(in.f10,glsthr,glat,zht,dtclst);

    // Compute the local exospheric temperature.
    // Add geomagnetic storm effect from input dTc value
    out.temp[0] = tsubl + in.dstdtc;
    const double tinf = tsubl + in.dstdtc + dtclst;

    // Equation 9
    const double tsubx = 444.3807e0 + 0.02385e0 * tinf -
                         392.8292e0 * std::exp(-0.0021357e0 * tinf);
    
    // Equation 11
    const double gsubx = 0.054285714e0 * (tsubx - 183e0);

    // The TC array will be an argument in the call to
    // XLOCAL, which evaluates Equation (10) or Equation (13)
    tc[0] = tsubx;
    tc[1] = gsubx;

    // A AND GSUBX/A OF Equation (13)
    tc[2] = (tinf - tsubx)/PIOV2;
    tc[3] = gsubx/tc[2];

    // Equation (5)
    const double z1 = 90e0;
    const double z2 = std::min(sat[2], 105e0);
    double al = std::log(z2 / z1);
    int n = static_cast<int>(al / r1) + 1;
    double zr = std::exp(al / static_cast<double>(n));
    const double ambar1 = xambar(z1);
    const double tloc1 = xlocal(z1, tc);
    double zend = z1;
    double sum2 = 0e0;
    double ain = ambar1 * xgrav(z1) / tloc1;

    double ambar2,tloc2,z,gravl;
    for (int i = 0; i < n; i++) {
      z = zend;
      zend = zr * z;
      const double dz = 0.25e0 * (zend - z);
      double sum1 = wt[0] * ain;
      for (int j = 1; j < 5; j++) {
        z = z + dz;
        ambar2 = xambar(z);
        tloc2 = xlocal(z, tc);
        const double gravl = xgrav(z);
        ain = ambar2 * gravl / tloc2;
        sum1 += wt[j] * ain;
      }
      sum2 += dz / sum1;
    }

    const double fact1 = 1000e0/rstar;
    out.rho =
        3.46e-6 * ambar2 * tloc1 * std::exp(-fact1 * sum2) / ambar1 / tloc2;

    // Equation (2)
    const double anm = avogad * out.rho;
    const double an  = anm/ambar2;

    // Equation (3)
    double fact2 = anm/28.960e0;
    aln[0] = std::log(frac[0]*fact2);
    aln[3] = std::log(frac[2]*fact2);
    aln[4] = std::log(frac[3]*fact2);

    // Equation (4)
    aln[1] = std::log(fact2 * (1e0 + frac[1]) - an);
    aln[2] = std::log(2e0 * (an - fact2));

  if (sat[2] <= 105e0) {
    out.temp[1] = tloc2;
    // Put in negligible hydrogen for use in DO-LOOP 13
    aln[5] = aln[4] - 25e0;
    // exit at bp 11
  } else {
    // (bp: 3) Equation (6)
    const double z3 = std::min(sat[2], 500e0);
    al = std::log(z3 / z);
    n = static_cast<int>(al / r2) + 1;
    zr = std::exp(al / (double)n);
    sum2 = 0e0;
    ain = gravl / tloc2;

    double tloc3;
    for (int i = 0; i < n; i++) {
      z = zend;
      zend = zr * z;
      const double dz = 0.25e0 * (zend - z);
      double sum1 = wt[0] * ain;
      for (int j = 1; j < 5; j++) {
        z = z + dz;
        tloc3 = xlocal(z, tc);
        gravl = xgrav(z);
        ain = gravl / tloc3;
        sum1 = sum1 + wt[j] * ain;
      }
      sum2 += dz * sum1
    }

    const double z4 = std::max(sat[2], 500e0);
    al = std::log(z4 / z);
    const double r = (sat[2] > 500e0) ? r3 : r2;
    n = static_cast<int>(al / r) + 1;
    zr = std::exp(al / (double)n);
    double sum3 = 0e0;

    double tloc4;
    for (int i=0; i<n; i++) {
      z = zend;
      zend = zr * z;
      const double dz = 0.25e0 * (zend - z);
      double sum1 = wt[0] * ain;
      for (int j=1; j<5; j++) {
        z = z + dz;
        tloc4 = xlocal(z,tc);
        gravl = xgrav(z);
        ain = gravl/tloc4;
        sum1 = sum1 + wt[j] * ain;
      }
      sum3 += dz * sum1;
    }

    double altr, hsigni, t500, hsign;
    if (sat[2] > 500e0) {
      t500 = tloc3;
      out.temp[1] = tloc4;
      altr = std::log(tloc4/tloc2);
      fact2 = fact1 * (sum2 + sum3);
      hsign = -1e0;
    } else {
      t500 = tloc4;
      out.temp[1] = tloc3;
      altr = std::log(tloc3/tloc2);
      fact2 = fact1 * sum2;
      hsign = 1e0;
    }

    for (int i = 0; i < 5; i++) {
      aln[i] -= (1e0 + alpha[i]) * altr - fact2 * amw[i];
    }
    // Equation (7) - Note that in CIRA72, AL10T5 = DLOG10(T500)
    const double al10t5 = std::log10(tinf);
    const double alnh5 = (5.5e0 * al10t5 - 39.40e0) * al10t5 + 73.13e0;
    aln[5] = al10 * (alnh5 + 6e0) +
             hsign * (std::log(tloc4 / tloc3) + fact1 * sum3 * amw[5]);
  }

  // bp: 11
  // Equation (24)  - J70 Seasonal-Latitudinal Variation
  const double trash = (amjd - 36204e0) / 365.2422e0;
  const double capphi = std::fmod(trash,1e0);

  const double dlrsl = 0.02e0 * (sat[2] - 90e0)
    * std::exp(-0.045e0 * (sat[2] - 90e0))
    * std::copysign(1e0,sat[1]) * std::sin(TWOPI * capphi + 1.72e0)
    * std::pow(std::sin(sat[1]), 2e0);

  // Equation (23) - Computes the semiannual variation
  double dlrsa = 0e0;
  if (z < 2000e0 ) {
    const double d1950 = amjd - 33281e0;
    tmoutd(d1950,iyr,yrday);
    // use new semiannual model
    double fzz,gtz;
    semian08(yrday,zht,in.f10b,in.s10b,in.xm10b,fzz,gtz,dlrsa);
    if (fzz < 0e0) dlrsa = 0e0;
  }

  // Sum the delta-log-rhos and apply to the number densities.
  // In CIRA72 the following equation contains an actual sum,
  // namely DLR = AL10 * (DLRGM + DLRSA + DLRSL)
  // However, for Jacchia 70, there is no DLRGM or DLRSA.
  const double dlr = al10 * (dlrsl + dlrsa);
  for (int i=0; i<6; i++) aln[i] += dlr;

  // Compute mass-density and mean-molecular-weight and
  // convert number density logs from natural to common.
  double sumn = 0e0, sumnm = 0e0;
  for (int i=0; i<6; i++) {
    const double an = std::exp(aln[i]);
    sumn = sumn + an;
    sumnm = sumnm + an*amw[i];
    al10n[i] = aln[i]/al10;
  }
  out.rho = sumnm/avogad;

  // Compute the high altitude exospheric density correction factor

  double fex = 1e0;
  if ((zht >= 1000e0) && (zht < 1500e0)) {
    const double zeta   = (zht - 1000e0) * 0.002e0;
    const double zeta2  =  zeta * zeta;
    const double zeta3  =  zeta * zeta2;
    const double f15c   = cht[0] + cht[1]*in.f10b + cht[2]*1500e0
      + cht[3]*in.f10b*1500e0;
    const double f15c_zeta = (cht[2] + cht[3]*in.f10b) * 500e0;
    const double fex2   = 3e0 * f15c - f15c_zeta - 3e0;
    const double fex3   = f15c_zeta - 2e0 * f15c + 2e0;
    fex    = 1e0 + fex2 * zeta2 + fex3 * zeta3;
  }
  if (zht >= 1500e0) {
    fex = cht[0] + cht[1]*in.f10b + cht[2]*zht + cht[3]*in.f10b*zht;
  }

  // Apply the exospheric density correction factor.
  out.rho *= fex;

  // all done
  return 0;
}
