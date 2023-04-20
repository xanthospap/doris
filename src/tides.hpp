#ifndef __IERS2010_VAR_TIDE_MODELS_HPP__
#define __IERS2010_VAR_TIDE_MODELS_HPP__

#include "iers2010/doodson.hpp"
#include "associated_legendre.hpp"
#include "base_error.hpp"
#include "cmat2d.hpp"
#include "datetime/dtcalendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "geodesy/units.hpp"
#include "stokes_coeffs.hpp"
#include <array>
#include <cstring>

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
  void clear_coefficients() noexcept;

#ifdef DEBUG
  void print_matrix_sizes() const noexcept;
#endif
}; // DoodsonOceanTideConstituent

int memmap_octide_coefficients(
    const char *fn, std::vector<dso::DoodsonOceanTideConstituent> &freqs,
    int max_degree, int max_order, int num_header_lines, double scale=1e0) noexcept;

class OceanTide {
private:
  std::vector<DoodsonOceanTideConstituent> doodsonFreqs;
  dso::StokesCoeffs dCS;
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> V; ///< workspace
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> W; ///< workspace
public:
  int max_degree() const noexcept {return dCS.max_degree();}
  int max_order() const noexcept {return dCS.max_order();}
  OceanTide(const std::vector<DoodsonOceanTideConstituent> &dfvec,
            double GMearth, double Rearth, int max_degree,
            int max_order) noexcept
      : doodsonFreqs(dfvec), dCS(max_degree, max_order, GMearth, Rearth),
        V(max_degree + 3, max_degree + 3), W(max_degree + 3, max_degree + 3) {}

  int operator()(const dso::TwoPartDate &mjdtt, int max_degree,
                 int max_order) noexcept;
  int acceleration(const dso::TwoPartDate &mjdtt,
                   const Eigen::Matrix<double, 3, 1> &rsat,
                   Eigen::Matrix<double, 3, 1> &acc, int max_degree = -1,
                   int max_order = -1) noexcept;
}; // OceanTide

class SolidEarthPoleTide {
  double m1,m2; // arcsec
 int operator()(const dso::TwoPartDate& mjdtt, double xp_sec, double yp_sec) noexcept {
   // secular pole (IERS2010, Sec7.2.4 Eq. 21)
   const double fyrs = mjdtt.as_fractional_years();
   const double xs = 55e0 + 1.677e0*(fyrs-2e3);    // [mas]
   const double ys = 320.5e0 + 3.460e0*(fyrs-2e3); // [mas]
   // "wobble" variables (IERS2010, Sec7.2.4 Eq. 25)
   m1 = xp_sec - xs*1e-3;    // [as]
   m2 = -(yp_sec - ys*1e-3); // [as]
   return 0;
 }
public:
 auto poleTide(const dso::TwoPartDate& mjdtt, double xp_sec, double yp_sec) noexcept {
   struct dCS21 { double dc21,ds21; };
   // compute "wobble" variables of date
   this->operator()(mjdtt,xp_sec,yp_sec);
   double dC21,dS21;
   // compute solid earth pole tide
   dC21 = -1.3331e-9*(m1+0.0115*m2);
   dS21 = -1.3331e-9*(m2-0.0115*m1);
   // ocean pole tide
   dC21 += -2.1778e-10*(m1-0.01724*m2);
   dS21 += -1.7232e-10*(m2-0.03365*m1);

   return dCS21{dC21,dS21};
 }
}; // SolidEarthPoleTide

class SolidEarthTide {
  static constexpr const int degree = 4;

private:
  const double GM_moon, GM_sun;
  dso::StokesCoeffs cs;
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
  int solid_earth_tide_step2_d(const dso::TwoPartDate &mjdtt,
                             double &dC20,
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
                   Eigen::Matrix<double, 3, 1> &acc,
                   Eigen::Matrix<double, 3, 3> &acc_gradient) noexcept;
  
  int acceleration(const dso::TwoPartDate &mjdtt,
                   const dso::TwoPartDate &mjdut1,
                   double xp_sec, double yp_sec,
                   const Eigen::Matrix<double, 3, 1> &rsat,
                   const Eigen::Matrix<double, 3, 1> &rmoon,
                   const Eigen::Matrix<double, 3, 1> &rsun,
                   Eigen::Matrix<double, 3, 1> &acc,
                   Eigen::Matrix<double, 3, 3> &acc_gradient) noexcept;
}; // SolidEarthTide

} // namespace dso

#endif
