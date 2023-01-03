#include "tides.hpp"

namespace {
/// @brief Third body (Sun or Moon) Solid earth tide geopotential coefficient
///        corrections, based on IERS 2010, Anelastic Earth
/// @param[in] Re Equatorial radius of the Earth [m]
/// @param[in] GM Gravitational constant of Earth
/// @param[in] r_tb ECEF, cartesian position vector of third body (sun or
///            Moon) [m]
/// @param[in] GM_tb Gravitational constant of Third body
/// @param[out] dC Array where the computed geopotential corrections ΔC are
///            added. Expected size is 12, in the order:
///             dC = C20,C21,C22,C30,C31,C32,C33,C40,C41,C42,C43,C44
///      [indexes] : 0   1   2   3   4   5   6   7   8   9   10  11
/// @param[out] dS Array where the computed geopotential corrections ΔS are
///            added. Expected size is 12, in the order:
///             dS = S20,S21,S22,0,S31,S32,S33,0,S41,S42,S43,S44
/// @return Always 0
int iers2010_solid_earth_tide_anelastic_tb(
    double Re, double GM, const Eigen::Matrix<double, 3, 1> &r_tb, double GM_tb,
    std::array<double, 12> &dC, std::array<double, 12> &dS) noexcept {
  // various variables
  const double rr = r_tb.norm(); // distance of third body to Earth's center
  const double sinphi = r_tb(2) / rr; // phi is φ (geocentric latitude)
  const double sinphi2 = sinphi * sinphi;
  const double cosphi = std::sqrt(1e0 - sinphi * sinphi);
  const double xlong = std::atan2(r_tb(1), r_tb(0)); // longitude of third body
  
  // compute associated Lagrange polynomials for n=2,3
  const double Pnm20 = std::sqrt(5e0) * 0.5e0 * (3e0 * sinphi2 - 1e0); // P20
  const double Pnm21 =
      std::sqrt(5.e0 / 3.e0) * 3.e0 * sinphi * std::sqrt(1.e0 - sinphi2); // P21
  const double Pnm22 = std::sqrt(5.e0 / 12.e0) * 3.e0 * (1.e0 - sinphi2); // P22
  const double Pnm30 = std::sqrt(7.e0) * 0.5e0 *
                       (5.e0 * sinphi2 * sinphi - 3.e0 * sinphi); // P30
  const double Pnm31 =
      std::sqrt(7.e0 / 6.e0) * 1.5e0 * (5.e0 * sinphi2 - 1.e0) * cosphi; // P31
  const double Pnm32 =
      std::sqrt(7.e0 / 60.e0) * 15.e0 * sinphi * cosphi * cosphi; // P32
  const double Pnm33 = std::sqrt(7.e0 / 360.e0) * 15.e0 * std::pow(cosphi, 3);

  // scale later on, watch only for the n=3 terms, which need an extra
  // scaling of (Re/r_j). Apply to all the 1/(2n+1) term
  constexpr const double fac2 = 1e0 / 5e0;
  const double fac3 = Re / rr / 7e0;

  // IERS 2010, Table 6.3: Nominal values of solid Earth
  // tide external potential Love numbers. Anelastic
  // Earth n  m  Re(k_nm)  Im(k_nm) k_nm^(+)
  // ----------------------------------
  // 2  0  0.30190    0.0
  // 2  1  0.29830   -0.00144
  // 2  2  0.30102   -0.00130
  // 3  0  0.09300
  // 3  1  0.09300
  // 3  2  0.09300
  // 3  3  0.09400
  const double cl = std::cos(xlong);
  const double sl = std::sin(xlong);
  const double c2l = std::cos(2 * xlong);
  const double s2l = std::sin(2 * xlong);
  const double c3l = std::cos(3 * xlong);
  const double s3l = std::sin(3 * xlong);

  // temporary storage
  std::array<double, 12> dCt = {0e0}, dSt = {0e0};

  // order n = 2
  dCt[0] = fac2 * Pnm20 * 0.30190e0;
  dCt[1] = fac2 * Pnm21 * (0.29830e0 * cl + (-0.00144e0) * sl);
  dSt[1] = fac2 * Pnm21 * (0.29830e0 * sl - (-0.00144e0) * cl);
  dCt[2] = fac2 * Pnm22 * (0.30102e0 * c2l + (-0.00130e0) * s2l);
  dSt[2] = fac2 * Pnm22 * (0.30102e0 * s2l - (-0.00130e0) * c2l);

  // order n = 3
  dCt[3] = fac3 * Pnm30 * 0.093e0;
  // dS30 = 0e0;
  dCt[4] = fac3 * Pnm31 * 0.093e0 * cl;
  dSt[4] = fac3 * Pnm31 * 0.093e0 * sl;
  dCt[5] = fac3 * Pnm32 * 0.093e0 * c2l;
  dSt[5] = fac3 * Pnm32 * 0.093e0 * s2l;
  dCt[6] = fac3 * Pnm33 * 0.094e0 * c3l;
  dSt[6] = fac3 * Pnm33 * 0.094e0 * s3l;

  // order n = 4
  dCt[7] = fac2 * Pnm20 * (-0.00089e0);
  // dS40 = 0e0;
  dCt[8] = fac2 * Pnm21 * (-0.00080e0) * cl;
  dSt[8] = fac2 * Pnm21 * (-0.00080e0) * sl;
  dCt[9] = fac2 * Pnm22 * (-0.00057e0) * c2l;
  dSt[9] = fac2 * Pnm22 * (-0.00057e0) * s2l;

  // scale and add to output arrays
  const double fac = (GM_tb / std::pow(rr, 3)) * (std::pow(Re, 3) / GM);
  std::transform(dCt.cbegin(), dCt.cend(), dC.cbegin(), dC.begin(),
                 [fac](double dct, double dc) { return dc + dct * fac; });
  std::transform(dSt.cbegin(), dSt.cend(), dS.cbegin(), dS.begin(),
                 [fac](double dst, double ds) { return ds + dst * fac; });

  return 0;
}
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
    const Eigen::Matrix<double, 3, 1> &rMoon,
    const Eigen::Matrix<double, 3, 1> &rSun, std::array<double, 12> &dC,
    std::array<double, 12> &dS) noexcept {

  // initialize corrections to zero
  std::fill(dC.begin(), dC.end(), 0e0);
  std::fill(dS.begin(), dS.end(), 0e0);

  // start with Sun geopotential corrections
  iers2010_solid_earth_tide_anelastic_tb(cs._Re, cs._GM, rSun, GM_sun, dC, dS);
  // add Moon
  iers2010_solid_earth_tide_anelastic_tb(cs._Re, cs._GM, rMoon, GM_moon, dC,
                                         dS);
  // all done for step 1
  return 0;
}
