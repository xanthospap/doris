#include "geodesy/geoconst.hpp"
#include "iers2010/iau.hpp"
#include "tides.hpp"
#include <cmath>
#include "geodesy/units.hpp"

namespace {
struct Step2EarthTideCoeffs {
  dso::DoodsonNumber d;
  ///< In-Phase Amp * 1e-12, Out-Of-Phase Amp * 1e-12
  double AIp, AOp;
}; // Step2EarthTideCoeffs

constexpr const Step2EarthTideCoeffs ST2_m20_D[]{
/*55565*/ {{0, 0, 0,  0,  1,  0},  1.660e+01, -6.700e+00},
/*55575*/ {{0, 0, 0,  0,  2,  0}, -1.000e-01,  1.000e-01},
/*56554*/ {{0, 0, 1,  0,  0, -1}, -1.200e+00,  8.000e-01},
/*57555*/ {{0, 0, 2,  0,  0,  0}, -5.500e+00,  4.300e+00},
/*57565*/ {{0, 0, 2,  0,  1,  0},  1.000e-01, -1.000e-01},
/*58554*/ {{0, 0, 3,  0,  0, -1}, -3.000e-01,  2.000e-01},
/*63655*/ {{0, 1,-2,  1,  0,  0}, -3.000e-01,  7.000e-01},
/*65445*/ {{0, 1, 0, -1, -1,  0},  1.000e-01, -2.000e-01},
/*65455*/ {{0, 1, 0, -1,  0,  0}, -1.200e+00,  3.700e+00},
/*65465*/ {{0, 1, 0, -1,  1,  0},  1.000e-01, -2.000e-01},
/*65655*/ {{0, 1, 0,  1,  0,  0},  1.000e-01, -2.000e-01},
/*73555*/ {{0, 2,-2,  0,  0,  0},  0.000e+00,  6.000e-01},
/*75355*/ {{0, 2, 0, -2,  0,  0},  0.000e+00,  3.000e-01},
/*75555*/ {{0, 2, 0,  0,  0,  0},  6.000e-01,  6.300e+00},
/*75565*/ {{0, 2, 0,  0,  1,  0},  2.000e-01,  2.600e+00},
/*75575*/ {{0, 2, 0,  0,  2,  0},  0.000e+00,  2.000e-01},
/*83655*/ {{0, 3,-2,  1,  0,  0},  1.000e-01,  2.000e-01},
/*85455*/ {{0, 3, 0, -1,  0,  0},  4.000e-01,  1.100e+00},
/*85465*/ {{0, 3, 0, -1,  1,  0},  2.000e-01,  5.000e-01},
/*93555*/ {{0, 4,-2,  0,  0,  0},  1.000e-01,  2.000e-01},
/*95355*/ {{0, 4, 0, -2,  0,  0},  1.000e-01,  1.000e-01}
};
constexpr const Step2EarthTideCoeffs ST2_m21_D[]{
/*125755*/ {{1, -3,  0,  2,  0,  0}, -10.00e-02, 00.00e+00},
/*127555*/ {{1, -3,  2,  0,  0,  0}, -10.00e-02, 00.00e+00},
/*135645*/ {{1, -2,  0,  1, -1,  0}, -10.00e-02, 00.00e+00},
/*135655*/ {{1, -2,  0,  1,  0,  0}, -70.00e-02, 10.00e-02},
/*137455*/ {{1, -2,  2, -1,  0,  0}, -10.00e-02, 00.00e+00},
/*145545*/ {{1, -1,  0,  0, -1,  0}, -01.30e+00, 10.00e-02},
/*145555*/ {{1, -1,  0,  0,  0,  0}, -06.80e+00, 60.00e-02},
/*147555*/ {{1, -1,  2,  0,  0,  0},  10.00e-02, 00.00e+00},
/*153655*/ {{1,  0, -2,  1,  0,  0},  10.00e-02, 00.00e+00},
/*155445*/ {{1,  0,  0, -1, -1,  0},  10.00e-02, 00.00e+00},
/*155455*/ {{1,  0,  0, -1,  0,  0},  40.00e-02, 00.00e+00},
/*155655*/ {{1,  0,  0,  1,  0,  0},  01.30e+00,-10.00e-02},
/*155665*/ {{1,  0,  0,  1,  1,  0},  30.00e-02, 00.00e+00},
/*157455*/ {{1,  0,  2, -1,  0,  0},  30.00e-02, 00.00e+00},
/*157465*/ {{1,  0,  2, -1,  1,  0},  10.00e-02, 00.00e+00},
/*162556*/ {{1,  1, -3,  0,  0,  1}, -01.90e+00, 10.00e-02},
/*163545*/ {{1,  1, -2,  0, -1,  0},  50.00e-02, 00.00e+00},
/*163555*/ {{1,  1, -2,  0,  0,  0}, -43.40e+00, 02.90e+00},
/*164554*/ {{1,  1, -1,  0,  0, -1},  60.00e-02, 00.00e+00},
/*164556*/ {{1,  1, -1,  0,  0,  1},  01.60e+00,-10.00e-02},
/*165345*/ {{1,  1,  0, -2, -1,  0},  10.00e-02, 00.00e+00},
/*165535*/ {{1,  1,  0,  0, -2,  0},  10.00e-02, 00.00e+00},
/*165545*/ {{1,  1,  0,  0, -1,  0}, -08.80e+00, 50.00e-02},
/*165555*/ {{1,  1,  0,  0,  0,  0},  04.71e+02,-30.20e+00},
/*165565*/ {{1,  1,  0,  0,  1,  0},  68.10e+00,-04.60e+00},
/*165575*/ {{1,  1,  0,  0,  2,  0}, -01.60e+00, 10.00e-02},
/*166455*/ {{1,  1,  1, -1,  0,  0},  10.00e-02, 00.00e+00},
/*166544*/ {{1,  1,  1,  0, -1, -1}, -10.00e-02, 00.00e+00},
/*166554*/ {{1,  1,  1,  0,  0, -1}, -20.60e+00,-30.00e-02},
/*166556*/ {{1,  1,  1,  0,  0,  1},  30.00e-02, 00.00e+00},
/*166564*/ {{1,  1,  1,  0,  1, -1}, -30.00e-02, 00.00e+00},
/*167355*/ {{1,  1,  2, -2,  0,  0}, -20.00e-02, 00.00e+00},
/*167365*/ {{1,  1,  2, -2,  1,  0}, -10.00e-02, 00.00e+00},
/*167555*/ {{1,  1,  2,  0,  0,  0}, -05.00e+00, 30.00e-02},
/*167565*/ {{1,  1,  2,  0,  1,  0},  20.00e-02, 00.00e+00},
/*168554*/ {{1,  1,  3,  0,  0, -1}, -20.00e-02, 00.00e+00},
/*173655*/ {{1,  2, -2,  1,  0,  0}, -50.00e-02, 00.00e+00},
/*173665*/ {{1,  2, -2,  1,  1,  0}, -10.00e-02, 00.00e+00},
/*175445*/ {{1,  2,  0, -1, -1,  0},  10.00e-02, 00.00e+00},
/*175455*/ {{1,  2,  0, -1,  0,  0}, -02.10e+00, 10.00e-02},
/*175465*/ {{1,  2,  0, -1,  1,  0}, -40.00e-02, 00.00e+00},
/*183555*/ {{1,  3, -2,  0,  0,  0}, -20.00e-02, 00.00e+00},
/*185355*/ {{1,  3,  0, -2,  0,  0}, -10.00e-02, 00.00e+00},
/*185555*/ {{1,  3,  0,  0,  0,  0}, -60.00e-02, 00.00e+00},
/*185565*/ {{1,  3,  0,  0,  1,  0}, -40.00e-02, 00.00e+00},
/*185575*/ {{1,  3,  0,  0,  2,  0}, -10.00e-02, 00.00e+00},
/*195455*/ {{1,  4,  0, -1,  0,  0}, -10.00e-02, 00.00e+00},
/*195465*/ {{1,  4,  0, -1,  1,  0}, -10.00e-02, 00.00e+00}};

constexpr const Step2EarthTideCoeffs ST2_m22_D[]{
    /*245655*/ {{2, -1, 0, 1, 0, 0}, -0.3e0, 0e0},
    /*255555*/ {{2, 0, 0, 0, 0, 0}, -1.2e0, 0e0}};

double compute_step2_m0_d(const double *const dargs) noexcept {
  constexpr const int szm20 = sizeof(ST2_m20_D) / sizeof(ST2_m20_D[0]);
  double dC20 = 0e0;
  for (int i = 0; i < szm20; i++) {
    const double theta = ST2_m20_D[i].d.phase(dargs);
    dC20 +=
        (ST2_m20_D[i].AIp * std::cos(theta) - ST2_m20_D[i].AOp * std::sin(theta));
  }
  return dC20 * 1e-12;
}
int compute_step2_m1_d(const double *const dargs, double &dC21,
                     double &dS21) noexcept {
  constexpr const int szm21 = sizeof(ST2_m21_D) / sizeof(ST2_m21_D[0]);
  // initial values for geopotential correction
  dC21 = 0e0;
  dS21 = 0e0;
  // iterate through Table 6.5a
  for (int i = 0; i < szm21; i++) {
    const double theta = ST2_m21_D[i].d.phase(dargs);
    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    dC21 += (ST2_m21_D[i].AIp * st + ST2_m21_D[i].AOp * ct);
    dS21 += (ST2_m21_D[i].AIp * ct - ST2_m21_D[i].AOp * st);
  }
  // mind the units!
  dC21 *= 1e-12;
  dS21 *= 1e-12;

  return 0;
}
int compute_step2_m2_d(const double *const dargs, double &dC22,
                     double &dS22) noexcept 
{
  constexpr const int szm22 = sizeof(ST2_m22_D) / sizeof(ST2_m22_D[0]);
  dC22 = 0e0;
  dS22 = 0e0;
  for (int i = 0; i < szm22; i++) {
    const double theta = ST2_m22_D[i].d.phase(dargs);
    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    dC22 += (ST2_m22_D[i].AIp * ct - ST2_m22_D[i].AOp * st);
    dS22 -= (ST2_m22_D[i].AIp * st - ST2_m22_D[i].AOp * ct);
  }
  // mind the units!
  dC22 *= 1e-12;
  dS22 *= 1e-12;

  return 0;
}
}// unnamed namespace
int dso::SolidEarthTide::solid_earth_tide_step2_d(const dso::TwoPartDate &mjdtt,
                                                double &dC20, double &dC21,
                                                double &dS21, double &dC22,
                                                double &dS22) const noexcept {

  // compute GMST **NOT** using IAU 2006/2000A [rad]
  const double gmst = dso::gmst_utc(mjdtt.tt2utc());

  // compute fundamental arguments (for given TT)
  const double fundarg[] = {
      iers2010::sofa::fal03(mjdtt),  // mean anomaly of moon, l
      iers2010::sofa::falp03(mjdtt), // mean anomaly of sun, l'
      iers2010::sofa::faf03(mjdtt),  // L - Ω, F
      iers2010::sofa::fad03(mjdtt),  // Mean Elongamjdttion of mjdtthe Moon from mjdtthe Sun, D
      iers2010::sofa::faom03(
          mjdtt) // Mean Longimjdttude of mjdtthe Ascending Node of mjdtthe Moon, Ω
  };

  // fundamental arguments to Doodson arguments
  double dargs[6];
  dso::fundarg2doodson(fundarg, gmst, dargs);

  // correction to C_20
  dC20 = compute_step2_m0_d(dargs);
  // correction to C_21 & S_21
  compute_step2_m1_d(dargs, dC21, dS21);
  // correction to C_22 & S_22
  compute_step2_m2_d(dargs, dC22, dS22);

  return 0;
}
