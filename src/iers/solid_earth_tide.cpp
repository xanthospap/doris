#include "tides.hpp"
#include "geodesy/geodesy.hpp"
#include "iers2010/iau.hpp"
#include <cmath>

namespace dso {
namespace
{
//    constexpr const double Am = 3.1274e-8; // [1/m] IERS, 6.8d

///< Nominal values of solid Earth tide external potential Love numbers.
///< IERS2010, Table 6.3
struct {
    int n,m;
    double knm, knm_Rplus, Reknm, Imknm, knm_Iplus;
} const LoveK[] = {{2, 0, 0.29525e0, -0.00087e0, 0.30190e0, 0e0, -0.00089e0},
             {2, 1, 0.29470e0, -0.00079e0, 0.29830e0, -0.00144e0, -0.00080e0},
             {2, 2, 0.29801e0, -0.00057e0, 0.30102e0, -0.00130e0, -0.00057e0},
             {3, 0, 0.093e0, 0e0, 0e0, 0e0, 0e0},
             {3,1,0.093e0, 0e0, 0e0, 0e0, 0e0},
             {3,2,0.093e0, 0e0, 0e0, 0e0, 0e0},
             {3,3,0.094e0, 0e0, 0e0, 0e0, 0e0}};

/// @brief < Table 6.5a from IERS2010, to compute Step 2 corrections for m=0
struct Step2TidesCoeffs {
    ///< Doodson Number is commented out
    ///< l,l’,F,D,Ω Multipliers for fundamental arguments
    int l,lp,F,D,Omega;
    ///< In-Phase Amp * 1e-12, Out-Of-Phase Amp * 1e-12
    double AIp, AOp;
}; // Step2Tides

/// @brief Table 6.5b from IERS2010. 
const Step2TidesCoeffs ST2_m20[] {
    {/*55565*/0,0,0,0,1,16.6e0,-6.7e0},
    {/*55575*/0,0,0,0,2,-0.1e0,0.1e0},
    {/*56554*/0,-1,0,0,0,-1.2e0,0.8e0},
    {/*57555*/0,0,-2,2,-2,-5.5e0,4.3e0},
    {/*57565*/0,0,-2,2,-1,0.1e0,-0.1e0},
    {/*58554*/0,-1,-2,2,-2,-0.3e0,0.2e0},
    {/*63655*/1,0,0,-2,0,-0.3e0,0.7e0},
    {/*65445*/-1,0,0,0,-1,0.1e0,-0.2e0},
    {/*65455*/-1,0,0,0,0,-1.2e0,3.7e0},
    {/*65465*/-1,0,0,0,1,0.1e0,-0.2e0},
    {/*65655*/1,0,-2,0,-2,0.1e0,-0.2e0},
    {/*73555*/0,0,0,-2,0,0.0e0,0.6e0},
    {/*75355*/-2,0,0,0,0,0.0e0,0.3e0},
    {/*75555*/0,0,-2,0,-2,0.6e0,6.3e0},
    {/*75565*/0,0,-2,0,-1,0.2e0,2.6e0},
    {/*75575*/0,0,-2,0,0,0.0e0,0.2e0},
    {/*83655*/1,0,-2,-2,-2,0.1e0,0.2e0},
    {/*85455*/-1,0,-2,0,-2,0.4e0,1.1e0},
    {/*85465*/-1,0,-2,0,-1,0.2e0,0.5e0},
    {/*93555*/0,0,-2,-2,-2,0.1e0,0.2e0},
    {/*95355*/-2,0,-2,0,-2,0.1e0,0.1e0},
};

/// @brief Table 6.5a from IERS2010. 
const Step2TidesCoeffs ST2_m21[] {
    {/*125755*/2,0,2,0,2,-0.1e0,0.0e0},
    {/*127555*/0,0,2,2,2,-0.1e0,0.0e0},
    {/*135645*/1,0,2,0,1,-0.1e0,0.0e0},
    {/*135655*/1,0,2,0,2,-0.7e0,0.1e0},
    {/*137455*/-1,0,2,2,2,-0.1e0,0.0e0},
    {/*145545*/0,0,2,0,1,-1.3e0,0.1e0},
    {/*145555*/0,0,2,0,2,-6.8e0,0.6e0},
    {/*147555*/0,0,0,2,0,0.1e0,0.0e0},
    {/*153655*/1,0,2,-2,2,0.1e0,0.0e0},
    {/*155445*/-1,0,2,0,1,0.1e0,0.0e0},
    {/*155455*/-1,0,2,0,2,0.4e0,0.0e0},
    {/*155655*/1,0,0,0,0,1.3e0,-0.1e0},
    {/*155665*/1,0,0,0,1,0.3e0,0.0e0},
    {/*157455*/-1,0,0,2,0,0.3e0,0.0e0},
    {/*157465*/-1,0,0,2,1,0.1e0,0.0e0},
    {/*162556*/0,1,2,-2,2,-1.9e0,0.1e0},
    {/*163545*/0,0,2,-2,1,0.5e0,0.0e0},
    {/*163555*/0,0,2,-2,2,-43.4e0,2.9e0},
    {/*164554*/0,-1,2,-2,2,0.6e0,0.0e0},
    {/*164556*/0,1,0,0,0,1.6e0,-0.1e0},
    {/*165345*/-2,0,2,0,1,0.1e0,0.0e0},
    {/*165535*/0,0,0,0,-2,0.1e0,0.0e0},
    {/*165545*/0,0,0,0,-1,-8.8e0,0.5e0},
    {/*165555*/0,0,0,0,0,470.9e0,-30.2e0},
    {/*165565*/0,0,0,0,1,68.1e0,-4.6e0},
    {/*165575*/0,0,0,0,2,-1.6e0,0.1e0},
    {/*166455*/-1,0,0,1,0,0.1e0,0.0e0},
    {/*166544*/0,-1,0,0,-1,-0.1e0,0.0e0},
    {/*166554*/0,-1,0,0,0,-20.6e0,-0.3e0},
    {/*166556*/0,1,-2,2,-2,0.3e0,0.0e0},
    {/*166564*/0,-1,0,0,1,-0.3e0,0.0e0},
    {/*167355*/-2,0,0,2,0,-0.2e0,0.0e0},
    {/*167365*/-2,0,0,2,1,-0.1e0,0.0e0},
    {/*167555*/0,0,-2,2,-2,-5.0e0,0.3e0},
    {/*167565*/0,0,-2,2,-1,0.2e0,0.0e0},
    {/*168554*/0,-1,-2,2,-2,-0.2e0,0.0e0},
    {/*173655*/1,0,0,-2,0,-0.5e0,0.0e0},
    {/*173665*/1,0,0,-2,1,-0.1e0,0.0e0},
    {/*175445*/-1,0,0,0,-1,0.1e0,0.0e0},
    {/*175455*/-1,0,0,0,0,-2.1e0,0.1e0},
    {/*175465*/-1,0,0,0,1,-0.4e0,0.0e0},
    {/*183555*/0,0,0,-2,0,-0.2e0,0.0e0},
    {/*185355*/-2,0,0,0,0,-0.1e0,0.0e0},
    {/*185555*/0,0,-2,0,-2,-0.6e0,0.0e0},
    {/*185565*/0,0,-2,0,-1,-0.4e0,0.0e0},
    {/*185575*/0,0,-2,0,0,-0.1e0,0.0e0},
    {/*195455*/-1,0,-2,0,-2,-0.1e0,0.0e0},
    {/*195465*/-1,0,-2,0,-1,-0.1e0,0.0e0}
};

/// @brief Compute Step 2 ΔC_(20) corrections, using Eq. 6.8a (IERS2010)
/// @param[in] fundarg Pointer to an array containing Fundamental Arguments,
///            in the order [l, l', F, F, Ω]
/// @return Step 2, ΔC_(20) correction
double compute_step2_m0(double gmst, const double const *fundarg) noexcept {
  constexpr const int szm20 = sizeof(ST2_m20) / sizeof(ST2_m20[0]);
    double dC20 = 0e0;
    // compute angle : θ = m*(θg + π) - Σ(N_j F_j), j=1,...5
    // Note that here m*(θ + π) is 0
    for (int i=0; i<szm20; i++) {
      const double theta = -(
          ST2_m20[i].l * fundarg[0] + ST2_m20[i].lp * fundarg[1] +
          ST2_m20[i].F * fundarg[2] + ST2_m20[i].D * fundarg[3] +
          ST2_m20[i].Omega * fundarg[4]);
      dC20 += (ST2_m20[i].AIp * std::cos(theta) -
              ST2_m20[i].AOp * std::sin(theta));
    }
    return dC20 * 1e-5;
}

/// @brief Compute Step 2, m=1 ΔC_(21) and ΔS_21 corrections, using Eq. 6.8b 
///        (IERS2010)
/// @param gmst Greenwich Mean Sidereal Time expressed in [rad]
/// @param[in] fundarg Pointer to an array containing Fundamental Arguments,
///                    in the order [l, l', F, F, Ω]
/// @return
int compute_step2_m1(double gmst, const double const *fundarg, double &dC21,
                     double &dS21) noexcept {
  constexpr const int szm21 = sizeof(ST2_m21) / sizeof(ST2_m21[0]);
  dC21 = 0e0;
  dS21 = 0e0;
  // iterate through Table 6.5a
  for (int i = 0; i < szm21; i++) {
    // compute angle : θ = m*(θg + π) - Σ(N_j F_j), j=1,...5
    const double theta =
        (gmst + dso::DPI / 2) -
        (ST2_m21[i].l * fundarg[0] + ST2_m21[i].lp * fundarg[1] +
         ST2_m21[i].F * fundarg[2] + ST2_m21[i].D * fundarg[3] +
         ST2_m21[i].Omega * fundarg[4]);
    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    dC21 += (ST2_m21[i].AIp * st + ST2_m21[i].AOp * ct);
    dS21 += (ST2_m21[i].AIp * ct - ST2_m21[i].AOp * st);
  }
  // mind the units!
  dC21 *= 1e-12;
  dS21 *= 1e-12;

  return 0;
}

/// @brief Compute Step 2, m=2 ΔC_(22) and ΔS_22 corrections, using Eq. 6.8b 
///        (IERS2010)
/// @param gmst Greenwich Mean Sidereal Time expressed in [rad]
/// @param[in] fundarg Pointer to an array containing Fundamental Arguments,
///                    in the order [l, l', F, F, Ω]
/// @return
int compute_step2_m2(double gmst, const double const *fundarg, double &dC22,
                     double &dS22) noexcept {
  dC22 = 0e0;
  dS22 = 0e0;
  
  // first for constituent N2 (245,655), see Table 6.5c
  // compute angle : θ = m*(θg + π) - Σ(N_j F_j), j=1,...5
  const double theta_n2 = (gmst + dso::DPI / 2) -
                          (1 * fundarg[0] + /*0 * fundarg[1]*/ +2 * fundarg[2] +
                           /*0 * fundarg[3]*/ +2 * fundarg[4]);
  // note that we only have corrections for the real part !
  dC22 += -0.3e0 * std::cos(theta_n2);
  dS22 -= -0.3e0 * std::sin(theta_n2);

  // second constituent M2
  const double theta_m2 =
      (gmst + dso::DPI / 2) - (/*0 * fundarg[0] + 0 * fundarg[1] +*/
                               2 * fundarg[2] + /*0 * fundarg[3] +*/
                               2 * fundarg[4]);
  // note that we only have corrections for the real part !
  dC22 += -1.2e0 * std::cos(theta_m2);
  dS22 -= -1.2e0 * std::sin(theta_m2);
  
  // mind the units!
  dC22 *= 1e-12;
  dS22 *= 1e-12;

  return 0;
}
} // namespace

struct SolidEarthTide {
    const double GM, GM_moon, GM_sun;
    const double Re;
    dso::datetime<dso::nanoseconds> t_tt, t_ut;
    dso::AssociatedLegendreFunctions p;

    /// @brief Compute the Step-1 effect of Solid Earth Tides, as in
    ///        IERS2010, Sec. 6.2.1 (Eq. 6.6 and Eq. 6.7)
    ///        Affects the ΔC_nm and ΔS_nm (correction) coefficients, for
    ///        (nm) = (20), (3,0), (4,0)
    ///               (21), (3,1), (4,1)
    ///               (22), (3,2), (4,2)
    ///                     (3,3)
    ///               -----|------|------
    ///                6.6   6.6    6.7      IERS2010 Equation
    /// @note It is expected that the Legendre polynomials passed in via the 
    ///       calling instance, are regularized
    /// @param[in] Rmoon Distance from geocenter to Moon [m]
    /// @param[in] Rsun  Distance from geocenter to Sun [m]
    /// @param[in] mlon  ECEF longitude (from Greenwich) of Moon [rad]
    /// @param[in] slon  ECEF longitude (from Greenwich) of Sun [rad]
    int solid_earth_tide_step1(double Rmoon, double Rsun, double mlon,
                               double slon) noexcept {
      const double RRm = Re/Rmoon;
      const double RRm3 = RRm *RRm *RRm;
      const double RRs = Re / Rsun;
      const double RRs3 = RRs * RRs * RRs;
      const double GMme = GM_moon / GM;
      const double GMse = GM_sun / GM;
      const double sml = std::sin(mlon);
      const double ssl = std::sin(slon);
      const double cml = std::cos(mlon);
      const double csl = std::cos(slon);
      const double s2ml = 2e0 * sml * cml;   //std::sin(2e0*mlon);
      const double s2sl = 2e0 * ssl * ssl;   //std::sin(2e0*slon);
      const double c2ml = 2e0*cml*cml - 1e0; //std::cos(2e0*mlon);
      const double c2sl = 2e0*csl*csl - 1e0; //std::cos(2e0*slon);

      // n = 2, m = 0
      double fac = (LoveK[0].knm) / (LoveK[0].n * 2 + 1);
      const double dc20_moon = fac * (GMme) * (RRm3) * p(2, 0);
      const double dc20_sun = fac *  (GMse) * (RRs3) * p(2, 0);

      // n = 2, m = 1
      fac = (LoveK[1].knm) / (LoveK[1].n * 2 + 1);
      const double dc21_moon = fac * (GMme) * (RRm3) * p(2,1) * cml;
      const double ds21_moon = fac * (GMme) * (RRm3) * p(2,1) * sml;
      const double dc21_sun =  fac * (GMse) * (RRs3) * p(2,1) * csl;
      const double ds21_sun =  fac * (GMse) * (RRs3) * p(2,1) * ssl;

      // n = 2, m = 2
      fac = (LoveK[2].knm) / (LoveK[2].n * 2 + 1);
      const double dc22_moon = fac * (GMme) * (RRm3) * p(2,2) * c2ml;
      const double ds22_moon = fac * (GMme) * (RRm3) * p(2,2) * s2ml;
      const double dc22_sun =  fac * (GMse) * (RRs3) * p(2,2) * c2sl;
      const double ds22_sun =  fac * (GMse) * (RRs3) * p(2,2) * s2sl;

      // n = 3, m = 0
      const double RRm4 = RRm3 * RRm;
      const double RRs4 = RRs3 * RRs;
      fac = (LoveK[3].knm) / (LoveK[3].n * 2 + 1);
      const double dc30_moon = fac * (GMme) * (RRm4) * p(3,0);
      const double dc30_sun  = fac * (GMse) * (RRs4) * p(3,0);
      
      // n = 3, m = 1
      fac = (LoveK[4].knm)/(LoveK[4].n*2+1);
      const double dc31_moon = fac * (GMme) * (RRm4) * p(3,1) * cml;
      const double ds31_moon = fac * (GMme) * (RRm4) * p(3,1) * sml;
      const double dc31_sun =  fac * (GMse) * (RRs4) * p(3,1) * csl;
      const double ds31_sun =  fac * (GMse) * (RRs4) * p(3,1) * ssl;

      // n = 3, m = 2
      fac = (LoveK[5].knm) / (LoveK[5].n * 2 + 1);
      const double dc32_moon = fac * (GMme) * (RRm4) * p(3,2) * c2ml;
      const double ds32_moon = fac * (GMme) * (RRm4) * p(3,2) * s2ml;
      const double dc32_sun =  fac * (GMse) * (RRs4) * p(3,2) * c2sl;
      const double ds32_sun =  fac * (GMse) * (RRs4) * p(3,2) * s2sl;
      
      // n = 3, m = 3
      fac = (LoveK[6].knm) / (LoveK[6].n * 2 + 1);
      const double dc33_moon = fac * (GMme) * (RRm4) * p(3,3) * std::cos(3e0*mlon);
      const double ds33_moon = fac * (GMme) * (RRm4) * p(3,3) * std::sin(3e0*mlon);
      const double dc33_sun =  fac * (GMse) * (RRs4) * p(3,3) * std::cos(3e0*slon);
      const double ds33_sun =  fac * (GMse) * (RRs4) * p(3,3) * std::sin(3e0*slon);
      
      // n = 4, m = 0
      fac = (LoveK[0].knm_Rplus) / 5;
      const double dc40_moon = fac * (GMme) * (RRm3) * p(2,0);
      const double dc40_sun  = fac * (GMse) * (RRs3) * p(2,0);
      
      // n = 4, m = 1
      fac = (LoveK[1].knm_Rplus) / 5e0;
      const double dc41_moon = fac * (GMme) * (RRm3) * p(2,1) * cml;
      const double ds41_moon = fac * (GMme) * (RRm3) * p(2,1) * sml;
      const double dc41_sun =  fac * (GMse) * (RRs3) * p(2,1) * csl;
      const double ds41_sun =  fac * (GMse) * (RRs3) * p(2,1) * ssl;

      // n = 4, m = 2
      fac = (LoveK[2].knm_Rplus) / 5e0;
      const double dc42_moon = fac * (GMme) * (RRm3) * p(2,2) * c2ml;
      const double ds42_moon = fac * (GMme) * (RRm3) * p(2,2) * s2ml;
      const double dc42_sun =  fac * (GMse) * (RRs3) * p(2,2) * c2sl;
      const double ds42_sun =  fac * (GMse) * (RRs3) * p(2,2) * s2sl;
      
      return 0;
    }

int solid_earth_tide_step2() const noexcept {

// compute GMST using IAU 2006/2000A
const double gmst = iers2010::sofa::gmst06(dso::mjd0_jd, t_ut.as_mjd(),
                                           dso::mjd0_jd, t_tt.as_mjd());

// compute Julian centuries since J2000.0 (TT)
const double t = t_tt.jcenturies_sinceJ2000();

// compute fundamental arguments
const double fundarg[] = {
    iers2010::sofa::fal03 (t), // mean anomaly of moon, l
    iers2010::sofa::falp03(t), // mean anomaly of sun, l'
    iers2010::sofa::faf03 (t), // L - Ω, F
    iers2010::sofa::fad03 (t), // Mean Elongation of the Moon from the Sun, D
    iers2010::sofa::faom03(t) // Mean Longitude of the Ascending Node of the Moon, Ω
};

const double dC20 = compute_step2_m0(gmst, fundarg);

double dC21, dS21;
compute_step2_m1(gmst, fundarg, dC21, dS21);

double dC22, dS22;
compute_step2_m2(gmst, fundarg, dC22, dS22);
}

int operator()(double mjd_tai) noexcept 
{
    p.compute(Moon_lat);
    const double Rmoon;

    p.compute(Sun_lat);
}
}; // SolidEarthTide

}// dso