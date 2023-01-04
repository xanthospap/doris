#ifndef __IERS2010_VAR_TIDE_MODELS_HPP__
#define __IERS2010_VAR_TIDE_MODELS_HPP__

#include "associated_legendre.hpp"
#include "base_error.hpp"
#include "cmat2d.hpp"
#include "datetime/dtcalendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "geodesy/units.hpp"
#include "harmonic_coeffs.hpp"
#include <array>
#include <cstring>

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
/// int doodson2intarray(const char *const str, int *arr) noexcept;

class DoodsonNumber {
private:
  int iar[6] = {0};

public:
  explicit DoodsonNumber(const char *str);
  explicit DoodsonNumber(const int *ar = nullptr) {
    if (ar)
      std::memcpy(iar, ar, sizeof(int) * 6);
  }
  char *str(char *buf) const noexcept;
  bool operator==(const DoodsonNumber &other) const noexcept {
    int diff = 0;
    for (int i = 0; i < 6; i++)
      diff += (iar[i] != other.iar[i]);
    return (diff == 0);
  }
  bool operator!=(const DoodsonNumber &other) const noexcept {
    return !(this->operator==(other));
  }
  int *multipliers(int *mults) const noexcept {
    mults[0] = iar[0];
    mults[1] = iar[1] - 5;
    mults[2] = iar[2] - 5;
    mults[3] = iar[3] - 5;
    mults[4] = iar[4] - 5;
    mults[5] = iar[5] - 5;
    return mults;
  }
};

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

class DoodsonOceanTideConstituent {
private:
  DoodsonNumber doodson;
  int maxl, maxm;
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *DelCpl{nullptr},
      *DelSpl{nullptr}, *DelCmi{nullptr}, *DelSmi{nullptr};
  void set_null() noexcept;
  void deallocate() noexcept;

public:
  DoodsonOceanTideConstituent() noexcept
      : doodson{}, maxl(0), maxm(0), DelCpl(nullptr), DelSpl(nullptr),
        DelCmi(nullptr), DelSmi(nullptr) {}
  DoodsonOceanTideConstituent(const dso::DoodsonNumber d, int max_degree,
                              int max_order);
  ~DoodsonOceanTideConstituent() noexcept;
  DoodsonOceanTideConstituent(const DoodsonOceanTideConstituent &) noexcept;
  DoodsonOceanTideConstituent(DoodsonOceanTideConstituent &&) noexcept;
  DoodsonOceanTideConstituent &
  operator=(const DoodsonOceanTideConstituent &) noexcept;
  DoodsonOceanTideConstituent &
  operator=(DoodsonOceanTideConstituent &&) noexcept;

  int max_degree() const noexcept { return maxl; }
  int max_order() const noexcept { return maxm; }
  const DoodsonNumber &doodson_number() const noexcept { return doodson; }
  /// @warning Your fault if *DelCpl is NULL!
  // const dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>& delCp()
  // const noexcept { return *DelCpl; } const
  // dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>& delSp() const
  // noexcept { return *DelSpl; } const
  // dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>& delCm() const
  // noexcept { return *DelCmi; } const
  // dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>& delSm() const
  // noexcept { return *DelSmi; }
  const double &delCp(int l, int m) const noexcept {
    return DelCpl->operator()(l, m);
  }
  const double &delSp(int l, int m) const noexcept {
    return DelSpl->operator()(l, m);
  }
  const double &delCm(int l, int m) const noexcept {
    return DelCmi->operator()(l, m);
  }
  const double &delSm(int l, int m) const noexcept {
    return DelSmi->operator()(l, m);
  }
  double &delCp(int l, int m) noexcept { return DelCpl->operator()(l, m); }
  double &delSp(int l, int m) noexcept { return DelSpl->operator()(l, m); }
  double &delCm(int l, int m) noexcept { return DelCmi->operator()(l, m); }
  double &delSm(int l, int m) noexcept { return DelSmi->operator()(l, m); }

  dso::iStatus resize(int maxDegree) noexcept;
#ifdef DEBUG
  void print_matrix_sizes() const noexcept;
#endif
}; // DoodsonOceanTideConstituent

/*int inspect_octide_coefficients(const char *fn,
 * std::vector<DoodsonOceanTideConstituent> &freqs) noexcept;*/
int memmap_octide_coefficients(
    const char *fn, std::vector<dso::DoodsonOceanTideConstituent> &freqs,
    int max_degree, int max_order, int num_header_lines) noexcept;
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
