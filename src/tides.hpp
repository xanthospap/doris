#ifndef __IERS2010_VAR_TIDE_MODELS_HPP__
#define __IERS2010_VAR_TIDE_MODELS_HPP__

#include "associated_legendre.hpp"
#include "datetime/dtcalendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "geodesy/units.hpp"
#include "harmonic_coeffs.hpp"
#include <array>

namespace dso {

/// @brief Resolve a Doodson number string to 6 integers.
/// These integers are actually the multipliers for the Doodson variables
/// [τ, s, h, p, N', ps], in this order.
/// @param[in] str A Doodson number as string. The string must be composed
///                of 3 ints, followed by a '.' or ',' character, followed
///                by 3 ints.
///                Ex. "055.555", "167,555"
/// @param[out] arr The resulting array of 6 integers
/// @return Anything other than 0 denotes an error
int doodson2intarray(const char *const str, int *arr) noexcept;

/// @brief Fundamental (Delaunay) arguments to Doodson variables.
/// All angles are in [rad] in the range [0,2π)
/// @param[in] fundarg Fundamental (Delaunay) arguments, in the order
///             [l, lp, f, d, Ω], see notes.
/// @param[out] doodson Corresponding Doodson variables, in the order
///             [τ, s, h, p, N', ps]
/// @note Explanation of symbols used:
///   * [0] l  : Mean anomaly of the Moon [rad]
///   * [1] lp : Mean anomaly of the Sun [rad]
///   * [2] f  : L - Ω [rad]
///   * [3] d  : Mean elongation of the Moon from the Sun [rad]
///   * [4] Ω  : Mean longitude of the ascending node of the Moon [rad]
///
///   * [0] WARNING TODO this element uses gmst
///   * [1] s  : Moon's mean longitude [rad]
///   * [2] h  : Sun's mean longitude [rad]
///   * [3] p  : Longitude of Moon's mean perigee
///   * [4] N' : Negative longitude of Moon's mean node
///   * [5] pl : Longitude of Sun's mean perigee
inline int fundarg2doodson(const double *const fundarg,
                           double *doodson) noexcept {
  doodson[1] = dso::anp(fundarg[2] + fundarg[4]);
  doodson[2] = dso::anp(fundarg[2] + fundarg[4] - fundarg[3]);
  doodson[3] = dso::anp(fundarg[2] + fundarg[4] - fundarg[0]);
  doodson[4] = dso::anp(-fundarg[4]);
  doodson[5] = dso::anp(fundarg[2] + fundarg[4] - fundarg[3] - fundarg[1]);
  return 0;
}

/// @brief Compute Pole Tides according to IERS 2010, Sec. 7.1.4
Eigen::Matrix<double, 3, 1>
pole_tide(double tfyears, double xp, double yp,
          const Eigen::Matrix<double, 3, 1> &r) noexcept;

/// @brief Compute Pole Tides according to IERS 2010, Sec. 7.1.4
int pole_tide(const dso::datetime<dso::nanoseconds> &t, double xp, double yp,
              const std::vector<Eigen::Matrix<double, 3, 1>> &sites,
              std::vector<Eigen::Matrix<double, 3, 1>> &Senu) noexcept;

/// @brief See Section 7.1.1.2 Permanent deformation, of IERS 2010
int permanent_tide(const std::vector<Eigen::Matrix<double, 3, 1>> &sites,
                   std::vector<Eigen::Matrix<double, 3, 1>> &Senu) noexcept;

class OceanTide {}; // OceanTide

class SolidEarthTide {
  static constexpr const int degree = 4;

private:
  const double GM_moon, GM_sun;
  dso::HarmonicCoeffs cs;
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> V; ///< workspace
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> W; ///< workspace
  // AssociatedLegendreFunctions PM, PS;

  int solid_earth_tide_step1(const Eigen::Matrix<double, 3, 1> &rmoon,
                             const Eigen::Matrix<double, 3, 1> &rsun,
                             std::array<double, 12> &dC,
                             std::array<double, 12> &dS) noexcept;

  int solid_earth_tide_step2(const dso::TwoPartDate &mjdtt,
                             const dso::TwoPartDate &mjdut1, double &dC20,
                             double &dC21, double &dS21, double &dC22,
                             double &dS22) const noexcept;

public:
  /// @brief Constructor
  /// @param GMearth Standard gravitational parameter μ=GM for the Earth
  /// [m^2/s^2]
  /// @param Rearth Equatorial radius of Earth [m]
  /// @param GMmoon Standard gravitational parameter μ=GM for the Moon
  ///        [m^2/s^2]
  /// @param GMsun Standard gravitational parameter μ=GM for the Moon
  ///        [m^2/s^2]
  SolidEarthTide(double GMearth, double Rearth, double GMmoon,
                 double GMsun) noexcept
      : GM_moon(GMmoon), GM_sun(GMsun), cs(degree, degree, GMearth, Rearth),
        V(degree + 3, degree + 3), W(degree + 3, degree + 3) /*, PM(degree),
         PS(degree)*/
  {}

  int operator()(const dso::TwoPartDate &mjdtt, const dso::TwoPartDate &mjdut1,
                 const Eigen::Matrix<double, 3, 1> &rmoon,
                 const Eigen::Matrix<double, 3, 1> &rsun,
                 std::array<double, 12> &dC,
                 std::array<double, 12> &dS) noexcept;

  int acceleration(const dso::TwoPartDate &mjdtt,
                   const dso::TwoPartDate &mjdut1,
                   const Eigen::Matrix<double, 3, 1> &rsat,
                   const Eigen::Matrix<double, 3, 1> &rmoon,
                   const Eigen::Matrix<double, 3, 1> &rsun,
                   Eigen::Matrix<double, 3, 1> &acc) noexcept;
}; // SolidEarthTide

} // namespace dso

#endif
