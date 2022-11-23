#include "tides.hpp"

namespace {
///< Nominal values of solid Earth tide external potential Love numbers.
///< IERS2010, Table 6.3
struct {
  int n, m;
  double knm;
} const LoveK[] = {
    {2, 0, 0.29525e0},
    {2, 1, 0.29470e0},
    {2, 2, 0.29801e0},
    {3, 0, 0.093e0},
    {3, 1, 0.093e0},
    {3, 2, 0.093e0},
    {3, 3, 0.094e0},
    {4, 0, -0.00087e0},   ///< k_20 ^(+), see Eq. 6.7
    {4, 1, -0.00079e0},   ///< k_21 ^(+)
    {4, 2, -0.00057e0}};  ///< k_22 ^(+)
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
  const double RRm = Re / Rmoon;
  const double RRm3 = RRm * RRm * RRm;
  const double RRs = Re / Rsun;
  const double RRs3 = RRs * RRs * RRs;
  const double GMme = GM_moon / GM;
  const double GMse = GM_sun / GM;
  const double sml = std::sin(mlon);
  const double ssl = std::sin(slon);
  const double cml = std::cos(mlon);
  const double csl = std::cos(slon);
  const double s2ml = 2e0 * sml * cml;       // std::sin(2e0*mlon);
  const double s2sl = 2e0 * ssl * ssl;       // std::sin(2e0*slon);
  const double c2ml = 2e0 * cml * cml - 1e0; // std::cos(2e0*mlon);
  const double c2sl = 2e0 * csl * csl - 1e0; // std::cos(2e0*slon);

  // n = 2, m = 0
  double fac = (LoveK[0].knm) / (LoveK[0].n * 2 + 1);
  const double dc20_moon = fac * (GMme) * (RRm3)*p(2, 0);
  const double dc20_sun = fac * (GMse) * (RRs3)*p(2, 0);

  // n = 2, m = 1
  fac = (LoveK[1].knm) / (LoveK[1].n * 2 + 1);
  const double dc21_moon = fac * (GMme) * (RRm3)*p(2, 1) * cml;
  const double ds21_moon = fac * (GMme) * (RRm3)*p(2, 1) * sml;
  const double dc21_sun = fac * (GMse) * (RRs3)*p(2, 1) * csl;
  const double ds21_sun = fac * (GMse) * (RRs3)*p(2, 1) * ssl;

  // n = 2, m = 2
  fac = (LoveK[2].knm) / (LoveK[2].n * 2 + 1);
  const double dc22_moon = fac * (GMme) * (RRm3)*p(2, 2) * c2ml;
  const double ds22_moon = fac * (GMme) * (RRm3)*p(2, 2) * s2ml;
  const double dc22_sun = fac * (GMse) * (RRs3)*p(2, 2) * c2sl;
  const double ds22_sun = fac * (GMse) * (RRs3)*p(2, 2) * s2sl;

  // n = 3, m = 0
  const double RRm4 = RRm3 * RRm;
  const double RRs4 = RRs3 * RRs;
  fac = (LoveK[3].knm) / (LoveK[3].n * 2 + 1);
  const double dc30_moon = fac * (GMme) * (RRm4)*p(3, 0);
  const double dc30_sun = fac * (GMse) * (RRs4)*p(3, 0);

  // n = 3, m = 1
  fac = (LoveK[4].knm) / (LoveK[4].n * 2 + 1);
  const double dc31_moon = fac * (GMme) * (RRm4)*p(3, 1) * cml;
  const double ds31_moon = fac * (GMme) * (RRm4)*p(3, 1) * sml;
  const double dc31_sun = fac * (GMse) * (RRs4)*p(3, 1) * csl;
  const double ds31_sun = fac * (GMse) * (RRs4)*p(3, 1) * ssl;

  // n = 3, m = 2
  fac = (LoveK[5].knm) / (LoveK[5].n * 2 + 1);
  const double dc32_moon = fac * (GMme) * (RRm4)*p(3, 2) * c2ml;
  const double ds32_moon = fac * (GMme) * (RRm4)*p(3, 2) * s2ml;
  const double dc32_sun = fac * (GMse) * (RRs4)*p(3, 2) * c2sl;
  const double ds32_sun = fac * (GMse) * (RRs4)*p(3, 2) * s2sl;

  // n = 3, m = 3
  fac = (LoveK[6].knm) / (LoveK[6].n * 2 + 1);
  const double dc33_moon = fac * (GMme) * (RRm4)*p(3, 3) * std::cos(3e0 * mlon);
  const double ds33_moon = fac * (GMme) * (RRm4)*p(3, 3) * std::sin(3e0 * mlon);
  const double dc33_sun = fac * (GMse) * (RRs4)*p(3, 3) * std::cos(3e0 * slon);
  const double ds33_sun = fac * (GMse) * (RRs4)*p(3, 3) * std::sin(3e0 * slon);

  // n = 4, m = 0
  fac = (LoveK[7].knm) / 5;
  const double dc40_moon = fac * (GMme) * (RRm3)*p(2, 0);
  const double dc40_sun = fac * (GMse) * (RRs3)*p(2, 0);

  // n = 4, m = 1
  fac = (LoveK[8].knm) / 5e0;
  const double dc41_moon = fac * (GMme) * (RRm3)*p(2, 1) * cml;
  const double ds41_moon = fac * (GMme) * (RRm3)*p(2, 1) * sml;
  const double dc41_sun = fac * (GMse) * (RRs3)*p(2, 1) * csl;
  const double ds41_sun = fac * (GMse) * (RRs3)*p(2, 1) * ssl;

  // n = 4, m = 2
  fac = (LoveK[9].knm) / 5e0;
  const double dc42_moon = fac * (GMme) * (RRm3)*p(2, 2) * c2ml;
  const double ds42_moon = fac * (GMme) * (RRm3)*p(2, 2) * s2ml;
  const double dc42_sun = fac * (GMse) * (RRs3)*p(2, 2) * c2sl;
  const double ds42_sun = fac * (GMse) * (RRs3)*p(2, 2) * s2sl;

  // assign corrections
  std::fill(std::begin(dC), std::end(dC), 0e0);
  dC[0] = dc20_moon + dc20_sun;
  dC[1] = dc21_moon + dc21_sun;
  dC[2] = dc22_moon + dc22_sun;
  dC[3] = dc30_moon + dc30_sun;
  dC[4] = dc31_moon + dc31_sun;
  dC[5] = dc32_moon + dc32_sun;
  dC[6] = dc33_moon + dc33_sun;
  dC[7] = dc40_moon + dc40_sun;
  dC[8] = dc41_moon + dc41_sun;
  dC[9] = dc42_moon + dc42_sun;

  std::fill(std::begin(dS), std::end(dS), 0e0);
  // dS[0] = 0e0;
  dS[1] = ds21_moon + ds21_sun;
  dS[2] = ds22_moon + ds22_sun;
  // dS[3] = 0e0
  dS[4] = ds31_moon + ds31_sun;
  dS[5] = ds32_moon + ds32_sun;
  dS[6] = ds33_moon + ds33_sun;
  // dS[7] =  0e0
  dS[8] = ds41_moon + ds41_sun;
  dS[9] = ds42_moon + ds42_sun;

  return 0;
}
