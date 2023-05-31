#include "geodesy/geoconst.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iau.hpp"
#include "tides.hpp"
#include <cmath>

namespace {

/* @brief Table 6.5a from IERS2010, to compute Step 2 corrections for m=0 */
struct Step2TidesCoeffs {
  /* Doodson Number is commented out
   * l,l’,F,D,Ω Multipliers for fundamental arguments
   */
  int l, lp, F, D, Omega;
  /* In-Phase Amp * 1e-12, Out-Of-Phase Amp * 1e-12 */
  double AIp, AOp;
}; /* Step2Tides */

/* @brief Table 6.5b from IERS2010.
 * Columns are: l l' F D Ω Amp(in-phase)*1e-12, Amp(out-of-phase)*1e-12
 * @warning Units in In- and Out-of- phase are 1e-12
 */
const Step2TidesCoeffs ST2_m20[]{
    {/*55565*/ 0, 0, 0, 0, 1, 16.6e0, -6.7e0},
    {/*55575*/ 0, 0, 0, 0, 2, -0.1e0, 0.1e0},
    {/*56554*/ 0, -1, 0, 0, 0, -1.2e0, 0.8e0},
    {/*57555*/ 0, 0, -2, 2, -2, -5.5e0, 4.3e0},
    {/*57565*/ 0, 0, -2, 2, -1, 0.1e0, -0.1e0},
    {/*58554*/ 0, -1, -2, 2, -2, -0.3e0, 0.2e0},
    {/*63655*/ 1, 0, 0, -2, 0, -0.3e0, 0.7e0},
    {/*65445*/ -1, 0, 0, 0, -1, 0.1e0, -0.2e0},
    {/*65455*/ -1, 0, 0, 0, 0, -1.2e0, 3.7e0},
    {/*65465*/ -1, 0, 0, 0, 1, 0.1e0, -0.2e0},
    {/*65655*/ 1, 0, -2, 0, -2, 0.1e0, -0.2e0},
    {/*73555*/ 0, 0, 0, -2, 0, 0.0e0, 0.6e0},
    {/*75355*/ -2, 0, 0, 0, 0, 0.0e0, 0.3e0},
    {/*75555*/ 0, 0, -2, 0, -2, 0.6e0, 6.3e0},
    {/*75565*/ 0, 0, -2, 0, -1, 0.2e0, 2.6e0},
    {/*75575*/ 0, 0, -2, 0, 0, 0.0e0, 0.2e0},
    {/*83655*/ 1, 0, -2, -2, -2, 0.1e0, 0.2e0},
    {/*85455*/ -1, 0, -2, 0, -2, 0.4e0, 1.1e0},
    {/*85465*/ -1, 0, -2, 0, -1, 0.2e0, 0.5e0},
    {/*93555*/ 0, 0, -2, -2, -2, 0.1e0, 0.2e0},
    {/*95355*/ -2, 0, -2, 0, -2, 0.1e0, 0.1e0},
};

/* @brief Table 6.5a from IERS2010.
 * Columns are: l l' F D Ω Amp(in-phase)*1e-12, Amp(out-of-phase)*1e-12
 * @warning Units in In- and Out-of- phase are 1e-12
 */
const Step2TidesCoeffs ST2_m21[]{{/*125755*/ 2, 0, 2, 0, 2, -0.1e0, 0.0e0},
                                 {/*127555*/ 0, 0, 2, 2, 2, -0.1e0, 0.0e0},
                                 {/*135645*/ 1, 0, 2, 0, 1, -0.1e0, 0.0e0},
                                 {/*135655*/ 1, 0, 2, 0, 2, -0.7e0, 0.1e0},
                                 {/*137455*/ -1, 0, 2, 2, 2, -0.1e0, 0.0e0},
                                 {/*145545*/ 0, 0, 2, 0, 1, -1.3e0, 0.1e0},
                                 {/*145555*/ 0, 0, 2, 0, 2, -6.8e0, 0.6e0},
                                 {/*147555*/ 0, 0, 0, 2, 0, 0.1e0, 0.0e0},
                                 {/*153655*/ 1, 0, 2, -2, 2, 0.1e0, 0.0e0},
                                 {/*155445*/ -1, 0, 2, 0, 1, 0.1e0, 0.0e0},
                                 {/*155455*/ -1, 0, 2, 0, 2, 0.4e0, 0.0e0},
                                 {/*155655*/ 1, 0, 0, 0, 0, 1.3e0, -0.1e0},
                                 {/*155665*/ 1, 0, 0, 0, 1, 0.3e0, 0.0e0},
                                 {/*157455*/ -1, 0, 0, 2, 0, 0.3e0, 0.0e0},
                                 {/*157465*/ -1, 0, 0, 2, 1, 0.1e0, 0.0e0},
                                 {/*162556*/ 0, 1, 2, -2, 2, -1.9e0, 0.1e0},
                                 {/*163545*/ 0, 0, 2, -2, 1, 0.5e0, 0.0e0},
                                 {/*163555*/ 0, 0, 2, -2, 2, -43.4e0, 2.9e0},
                                 {/*164554*/ 0, -1, 2, -2, 2, 0.6e0, 0.0e0},
                                 {/*164556*/ 0, 1, 0, 0, 0, 1.6e0, -0.1e0},
                                 {/*165345*/ -2, 0, 2, 0, 1, 0.1e0, 0.0e0},
                                 {/*165535*/ 0, 0, 0, 0, -2, 0.1e0, 0.0e0},
                                 {/*165545*/ 0, 0, 0, 0, -1, -8.8e0, 0.5e0},
                                 {/*165555*/ 0, 0, 0, 0, 0, 470.9e0, -30.2e0},
                                 {/*165565*/ 0, 0, 0, 0, 1, 68.1e0, -4.6e0},
                                 {/*165575*/ 0, 0, 0, 0, 2, -1.6e0, 0.1e0},
                                 {/*166455*/ -1, 0, 0, 1, 0, 0.1e0, 0.0e0},
                                 {/*166544*/ 0, -1, 0, 0, -1, -0.1e0, 0.0e0},
                                 {/*166554*/ 0, -1, 0, 0, 0, -20.6e0, -0.3e0},
                                 {/*166556*/ 0, 1, -2, 2, -2, 0.3e0, 0.0e0},
                                 {/*166564*/ 0, -1, 0, 0, 1, -0.3e0, 0.0e0},
                                 {/*167355*/ -2, 0, 0, 2, 0, -0.2e0, 0.0e0},
                                 {/*167365*/ -2, 0, 0, 2, 1, -0.1e0, 0.0e0},
                                 {/*167555*/ 0, 0, -2, 2, -2, -5.0e0, 0.3e0},
                                 {/*167565*/ 0, 0, -2, 2, -1, 0.2e0, 0.0e0},
                                 {/*168554*/ 0, -1, -2, 2, -2, -0.2e0, 0.0e0},
                                 {/*173655*/ 1, 0, 0, -2, 0, -0.5e0, 0.0e0},
                                 {/*173665*/ 1, 0, 0, -2, 1, -0.1e0, 0.0e0},
                                 {/*175445*/ -1, 0, 0, 0, -1, 0.1e0, 0.0e0},
                                 {/*175455*/ -1, 0, 0, 0, 0, -2.1e0, 0.1e0},
                                 {/*175465*/ -1, 0, 0, 0, 1, -0.4e0, 0.0e0},
                                 {/*183555*/ 0, 0, 0, -2, 0, -0.2e0, 0.0e0},
                                 {/*185355*/ -2, 0, 0, 0, 0, -0.1e0, 0.0e0},
                                 {/*185555*/ 0, 0, -2, 0, -2, -0.6e0, 0.0e0},
                                 {/*185565*/ 0, 0, -2, 0, -1, -0.4e0, 0.0e0},
                                 {/*185575*/ 0, 0, -2, 0, 0, -0.1e0, 0.0e0},
                                 {/*195455*/ -1, 0, -2, 0, -2, -0.1e0, 0.0e0},
                                 {/*195465*/ -1, 0, -2, 0, -1, -0.1e0, 0.0e0}};

/* @brief Compute Step 2 ΔC_(20) corrections, using Eq. 6.8a (IERS2010)
 * @param[in] fundarg Pointer to an array containing Fundamental Arguments,
 *            in the order [l, l', F, F, Ω]
 * @return Step 2, ΔC_(20) correction
 */
double compute_step2_m0([[maybe_unused]] double gmst,
                        const double *const fundarg) noexcept {
  constexpr const int szm20 = sizeof(ST2_m20) / sizeof(ST2_m20[0]);
  double dC20 = 0e0;
  /* compute angle : θ = m*(θg + π) - Σ(N_j F_j), j=1,...5
   * Note that here m*(θ + π) is 0
   */
  for (int i = 0; i < szm20; i++) {
    const double theta =
        dso::anp(-(ST2_m20[i].l * fundarg[0] + ST2_m20[i].lp * fundarg[1] +
                   ST2_m20[i].F * fundarg[2] + ST2_m20[i].D * fundarg[3] +
                   ST2_m20[i].Omega * fundarg[4]));
    dC20 +=
        (ST2_m20[i].AIp * std::cos(theta) - ST2_m20[i].AOp * std::sin(theta));
  }
  return dC20 * 1e-12;
}

/* @brief Compute Step 2, m=1 ΔC_(21) and ΔS_21 corrections, using Eq. 6.8b
 *        (IERS2010)
 * @param gmst Greenwich Mean Sidereal Time expressed in [rad]
 * @param[in] fundarg Pointer to an array containing Fundamental Arguments,
 *                    in the order [l, l', F, F, Ω]
 * @return
 */
int compute_step2_m1(double gmst, const double *const fundarg, double &dC21,
                     double &dS21) noexcept {
  constexpr const int szm21 = sizeof(ST2_m21) / sizeof(ST2_m21[0]);
  /* θ_g + π in [rad] */
  const double g = dso::anp(gmst + dso::DPI);
  /* initial values for geopotential correction */
  dC21 = 0e0;
  dS21 = 0e0;

  /* iterate through Table 6.5a */
  for (int i = 0; i < szm21; i++) {
    /* compute angle : θ = m * (θg + π) - Σ(N_j F_j), j=1,...5 and here m=1 */
    const double theta =
        dso::anp(g - (ST2_m21[i].l * fundarg[0] + ST2_m21[i].lp * fundarg[1] +
                      ST2_m21[i].F * fundarg[2] + ST2_m21[i].D * fundarg[3] +
                      ST2_m21[i].Omega * fundarg[4]));
    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    dC21 += (ST2_m21[i].AIp * st + ST2_m21[i].AOp * ct);
    dS21 += (ST2_m21[i].AIp * ct - ST2_m21[i].AOp * st);
  }

  /* mind the units! */
  dC21 *= 1e-12;
  dS21 *= 1e-12;

  return 0;
}

/* @brief Compute Step 2, m=2 ΔC_(22) and ΔS_22 corrections, using Eq. 6.8b
 *        (IERS2010)
 * @param gmst Greenwich Mean Sidereal Time expressed in [rad]
 * @param[in] fundarg Pointer to an array containing Fundamental Arguments,
 *                    in the order [l, l', F, F, Ω]
 * @return
 */
int compute_step2_m2(double gmst, const double *const fundarg, double &dC22,
                     double &dS22) noexcept {
  /* θ_g + π in [rad] */
  const double g = dso::anp(gmst + dso::DPI);
  dC22 = 0e0;
  dS22 = 0e0;

  /* first for constituent N2 (245,655), see Table 6.5c
   * compute angle : θ = m*(θg + π) - Σ(N_j F_j), j=1,...5, here m=2
   */
  const double theta_n2 =
      dso::anp(2e0 * g - (1 * fundarg[0] + /*0 * fundarg[1]*/ +2 * fundarg[2] +
                          /*0 * fundarg[3]*/ +2 * fundarg[4]));

  /* note that we only have corrections for the real part ! */
  dC22 += -0.3e0 * std::cos(theta_n2);
  dS22 -= -0.3e0 * std::sin(theta_n2);

  /* second constituent M2 */
  const double theta_m2 =
      dso::anp(2e0 * g - (                 /*0 * fundarg[0] + 0 * fundarg[1] +*/
                          2 * fundarg[2] + /*0 * fundarg[3] +*/
                          2 * fundarg[4]));

  /* note that we only have corrections for the real part ! */
  dC22 += -1.2e0 * std::cos(theta_m2);
  dS22 -= -1.2e0 * std::sin(theta_m2);

  /* mind the units! */
  dC22 *= 1e-12;
  dS22 *= 1e-12;

  return 0;
}
} /* unnamed namespace */

int dso::SolidEarthTide::solid_earth_tide_step2(const dso::TwoPartDate &mjdtt,
                                                const dso::TwoPartDate &mjdut1,
                                                double &dC20, double &dC21,
                                                double &dS21, double &dC22,
                                                double &dS22) const noexcept {

  /* compute GMST using IAU 2006/2000A [rad] */
  const double gmst = iers2010::sofa::gmst06(mjdut1, mjdtt);

  /* compute fundamental arguments (for given TT) */
  const double fundarg[] = {iers2010::fal03(mjdtt), iers2010::falp03(mjdtt),
                            iers2010::faf03(mjdtt), iers2010::fad03(mjdtt),
                            iers2010::faom03(mjdtt)};

  /* correction to C_20 */
  dC20 = compute_step2_m0(gmst, fundarg);
  /* correction to C_21 & S_21 */
  compute_step2_m1(gmst, fundarg, dC21, dS21);
  /* correction to C_22 & S_22 */
  compute_step2_m2(gmst, fundarg, dC22, dS22);

  return 0;
}
