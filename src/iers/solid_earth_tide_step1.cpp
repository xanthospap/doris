#include "tides.hpp"

namespace {
///< Nominal values of solid Earth tide external potential Love numbers.
///< IERS2010, Table 6.3
struct {
  int n, m;
  double knm;
} const LoveK[] = {
    {2, 0, 0.29525e0},  {2, 1, 0.29470e0},  {2, 2, 0.29801e0},
    {3, 0, 0.093e0},    {3, 1, 0.093e0},    {3, 2, 0.093e0},
    {3, 3, 0.094e0},    {4, 0, -0.00087e0}, ///< k_20 ^(+), see Eq. 6.7
    {4, 1, -0.00079e0},                     ///< k_21 ^(+)
    {4, 2, -0.00057e0}};                    ///< k_22 ^(+)
} // unnamed namespace

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
/// @param[out] dC   Normalized corrections to C coefficients, in the order:
///             dC = C20,C21,C22,C30,C31,C32,C33,C40,C41,C42,0,0
/// @param[out] dS   Normalized corrections to S coefficients, in the order:
///             dS = 0,S21,S22,0,S31,S32,S33,0,S41,S42,0,0
int dso::SolidEarthTide::solid_earth_tide_step1(
    double Rmoon, double Rsun, double mlon, double slon,
    std::array<double, 12> &dC, std::array<double, 12> &dS) noexcept {
  const double Re = cs._Re;
  const double GM = cs._GM;

  const double RRm = Re / Rmoon;
  const double RRs = Re / Rsun;
  // do not use this factor, scale results later on
  const double RRm3 = 1e0; // RRm * RRm * RRm;
  // do not use this factor, scale results later on
  const double RRs3 = 1e0; // RRs * RRs * RRs;
  // do not use this factor, scale results later on
  const double GMme = 1e0; // GM_moon / GM;
  // do not use this factor, scale results later on
  const double GMse = 1e0; // GM_sun / GM;
  const double sml = std::sin(mlon);
  const double ssl = std::sin(slon);
  const double cml = std::cos(mlon);
  const double csl = std::cos(slon);
  const double s2ml = 2e0 * sml * cml;       // std::sin(2e0*mlon);
  const double s2sl = 2e0 * ssl * ssl;       // std::sin(2e0*slon);
  const double c2ml = 2e0 * cml * cml - 1e0; // std::cos(2e0*mlon);
  const double c2sl = 2e0 * csl * csl - 1e0; // std::cos(2e0*slon);

  // n = 2, m = 0, fac = k_20 / 2*2+1
  double fac = (LoveK[0].knm) / (LoveK[0].n * 2e0 + 1);
  const double dc20_moon = fac * GMme * RRm3 * PM(2, 0);
  //const double dc20_sun = fac * GMse * RRs3 * PS(2, 0);
  const double dc20_sun =
      PS(2, 0) * 0.30190e0 * cml / 5e0;

  // n = 2, m = 1, fac = k_21 / 2*2+1
  fac = (LoveK[1].knm) / (LoveK[1].n * 2 + 1);
  const double dc21_moon = fac * GMme * RRm3 * PM(2, 1) * cml;
  const double ds21_moon = fac * GMme * RRm3 * PM(2, 1) * sml;
  const double dc21_sun = fac * GMse * RRs3 * PS(2, 1) * csl;
  const double ds21_sun = fac * GMse * RRs3 * PS(2, 1) * ssl;

  // n = 2, m = 2, fac = k_22 / 2*2+1
  fac = (LoveK[2].knm) / (LoveK[2].n * 2 + 1);
  const double dc22_moon = fac * GMme * RRm3 * PM(2, 2) * c2ml;
  const double ds22_moon = fac * GMme * RRm3 * PM(2, 2) * s2ml;
  const double dc22_sun = fac * GMse * RRs3 * PS(2, 2) * c2sl;
  const double ds22_sun = fac * GMse * RRs3 * PS(2, 2) * s2sl;

  // n = 3, m = 0, fac = k_30 / 2*3+1
  const double RRm4 = RRm3 * RRm;
  const double RRs4 = RRs3 * RRs;
  fac = (LoveK[3].knm) / (LoveK[3].n * 2 + 1);
  const double dc30_moon = fac * GMme * RRm4 * PM(3, 0);
  const double dc30_sun = fac * GMse * RRs4 * PS(3, 0);

  // n = 3, m = 1, fac = k_31 / 2*3+1
  fac = (LoveK[4].knm) / (LoveK[4].n * 2 + 1);
  const double dc31_moon = fac * GMme * RRm4 * PM(3, 1) * cml;
  const double ds31_moon = fac * GMme * RRm4 * PM(3, 1) * sml;
  const double dc31_sun = fac * GMse * RRs4 * PS(3, 1) * csl;
  const double ds31_sun = fac * GMse * RRs4 * PS(3, 1) * ssl;

  // n = 3, m = 2, fac = k_32 / 2*3+1
  fac = (LoveK[5].knm) / (LoveK[5].n * 2 + 1);
  const double dc32_moon = fac * GMme * RRm4 * PM(3, 2) * c2ml;
  const double ds32_moon = fac * GMme * RRm4 * PM(3, 2) * s2ml;
  const double dc32_sun = fac * GMse * RRs4 * PS(3, 2) * c2sl;
  const double ds32_sun = fac * GMse * RRs4 * PS(3, 2) * s2sl;

  // n = 3, m = 3, fac = k_33 / 2*3+1
  fac = (LoveK[6].knm) / (LoveK[6].n * 2 + 1);
  const double dc33_moon = fac * GMme * RRm4 * PM(3, 3) * std::cos(3e0 * mlon);
  const double ds33_moon = fac * GMme * RRm4 * PM(3, 3) * std::sin(3e0 * mlon);
  const double dc33_sun = fac * GMse * RRs4 * PS(3, 3) * std::cos(3e0 * slon);
  const double ds33_sun = fac * GMse * RRs4 * PS(3, 3) * std::sin(3e0 * slon);

  // n = 4, m = 0, fac = k_40 / 5
  fac = (LoveK[7].knm) / 5e0;
  const double dc40_moon = fac * GMme * RRm3 * PM(2, 0);
  const double dc40_sun = fac * GMse * RRs3 * PS(2, 0);

  // n = 4, m = 1, fac = k_41 / 5
  fac = (LoveK[8].knm) / 5e0;
  const double dc41_moon = fac * GMme * RRm3 * PM(2, 1) * cml;
  const double ds41_moon = fac * GMme * RRm3 * PM(2, 1) * sml;
  const double dc41_sun = fac * GMse * RRs3 * PS(2, 1) * csl;
  const double ds41_sun = fac * GMse * RRs3 * PS(2, 1) * ssl;

  // n = 4, m = 2, fac = k_42 / 5
  fac = (LoveK[9].knm) / 5e0;
  const double dc42_moon = fac * GMme * RRm3 * PM(2, 2) * c2ml;
  const double ds42_moon = fac * GMme * RRm3 * PM(2, 2) * s2ml;
  const double dc42_sun = fac * GMse * RRs3 * PS(2, 2) * c2sl;
  const double ds42_sun = fac * GMse * RRs3 * PS(2, 2) * s2sl;

  // scale
  const double ReRm3 = RRm * RRm * RRm;
  const double ReRs3 = RRs * RRs * RRs;
  const double GmMGmE = GM_moon / GM;
  const double GmSGmE = GM_sun / GM;

  // assign corrections
  std::fill(std::begin(dC), std::end(dC), 0e0);
  dC[0] = dc20_moon * (GmMGmE * ReRm3) + dc20_sun * (GmSGmE * ReRs3);
  dC[1] = dc21_moon * (GmMGmE * ReRm3) + dc21_sun * (GmSGmE * ReRs3);
  dC[2] = dc22_moon * (GmMGmE * ReRm3) + dc22_sun * (GmSGmE * ReRs3);
  dC[3] = dc30_moon * (GmMGmE * ReRm3) + dc30_sun * (GmSGmE * ReRs3);
  dC[4] = dc31_moon * (GmMGmE * ReRm3) + dc31_sun * (GmSGmE * ReRs3);
  dC[5] = dc32_moon * (GmMGmE * ReRm3) + dc32_sun * (GmSGmE * ReRs3);
  dC[6] = dc33_moon * (GmMGmE * ReRm3) + dc33_sun * (GmSGmE * ReRs3);
  dC[7] = dc40_moon * (GmMGmE * ReRm3) + dc40_sun * (GmSGmE * ReRs3);
  dC[8] = dc41_moon * (GmMGmE * ReRm3) + dc41_sun * (GmSGmE * ReRs3);
  dC[9] = dc42_moon * (GmMGmE * ReRm3) + dc42_sun * (GmSGmE * ReRs3);
  printf("\tC(%d,%d) = %.15e (P(2,0)*cml=%.12e)\n", 2,0,dc20_sun * (GmSGmE * ReRs3), PS(2,0)*cml);
  printf("\tC(%d,%d) = %.15e\n", 2,1,dc21_sun * (GmSGmE * ReRs3));
  printf("\tC(%d,%d) = %.15e\n", 2,2,dc22_sun * (GmSGmE * ReRs3));
  printf("\tC(%d,%d) = %.15e\n", 3,0,dc30_sun * (GmSGmE * ReRs3));
  printf("\tC(%d,%d) = %.15e\n", 3,1,dc31_sun * (GmSGmE * ReRs3));
  printf("\tC(%d,%d) = %.15e\n", 3,2,dc32_sun * (GmSGmE * ReRs3));
  printf("\tC(%d,%d) = %.15e\n", 3,3,dc33_sun * (GmSGmE * ReRs3));
  printf("\tC(%d,%d) = %.15e\n", 4,0,dc40_sun * (GmSGmE * ReRs3));
  printf("\tC(%d,%d) = %.15e\n", 4,1,dc41_sun * (GmSGmE * ReRs3));
  printf("\tC(%d,%d) = %.15e\n", 4,2,dc42_sun * (GmSGmE * ReRs3));

  std::fill(std::begin(dS), std::end(dS), 0e0);
  // dS[0] = 0e0;
  dS[1] = ds21_moon * (GmMGmE * ReRm3) + ds21_sun * (GmSGmE * ReRs3);
  dS[2] = ds22_moon * (GmMGmE * ReRm3) + ds22_sun * (GmSGmE * ReRs3);
  // dS[3] = 0e0
  dS[4] = ds31_moon * (GmMGmE * ReRm3) + ds31_sun * (GmSGmE * ReRs3);
  dS[5] = ds32_moon * (GmMGmE * ReRm3) + ds32_sun * (GmSGmE * ReRs3);
  dS[6] = ds33_moon * (GmMGmE * ReRm3) + ds33_sun * (GmSGmE * ReRs3);
  // dS[7] =  0e0
  dS[8] = ds41_moon * (GmMGmE * ReRm3) + ds41_sun * (GmSGmE * ReRs3);
  dS[9] = ds42_moon * (GmMGmE * ReRm3) + ds42_sun * (GmSGmE * ReRs3);

  return 0;
}
