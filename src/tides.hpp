#ifndef __IERS2010_VAR_TIDE_MODELS_HPP__
#define __IERS2010_VAR_TIDE_MODELS_HPP__

#include "associated_legendre.hpp"
#include "datetime/dtcalendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "harmonic_coeffs.hpp"
#include <array>

namespace dso {

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
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> V;
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> W;
  AssociatedLegendreFunctions PM, PS;

  int solid_earth_tide_step1(double Rmoon, double Rsun, double mlon,
                             double slon, std::array<double, 12> &dC,
                             std::array<double, 12> &dS) noexcept;

  int solid_earth_tide_step2(const dso::datetime<dso::nanoseconds> &t_tt,
                             double ut1_mjd, double &dC20, double &dC21,
                             double &dS21, double &dC22,
                             double &dS22) const noexcept;

public:
  /// @brief Constructor
  /// @param GMearth Standard gravitational parameter μ=GM for the Earth
  /// [m^2/s^2]
  /// @param Rearth Equatorial radius of Earth [m]
  /// @param GMmoon Standard gravitational parameter μ=GM for the Moon [m^2/s^2]
  /// @param GMsun Standard gravitational parameter μ=GM for the Moon [m^2/s^2]
  SolidEarthTide(double GMearth, double Rearth, double GMmoon,
                 double GMsun) noexcept
      : GM_moon(GMmoon), GM_sun(GMsun), cs(degree, degree, GMearth, Rearth),
        V(degree + 3, degree + 3), W(degree + 3, degree + 3), PM(degree),
        PS(degree) {}

  int operator()(/*dso::datetime<dso::nanoseconds> &t_tt,
                 double ut1_mjd,*/
                 const Eigen::Matrix<double, 3, 1> &rmoon,
                 const Eigen::Matrix<double, 3, 1> &rsun,
                 std::array<double, 12> &dC,
                 std::array<double, 12> &dS) noexcept;

  int acceleration(const Eigen::Matrix<double, 3, 1> &rsat,
                   const Eigen::Matrix<double, 3, 1> &rmoon,
                   const Eigen::Matrix<double, 3, 1> &rsun,
                   Eigen::Matrix<double, 3, 1> &acc) noexcept;
}; // SolidEarthTide

} // namespace dso

#endif
