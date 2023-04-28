#ifndef __IERS2010_VAR_TIDE_MODELS_HPP__
#define __IERS2010_VAR_TIDE_MODELS_HPP__

#include "associated_legendre.hpp"
#include "base_error.hpp"
#include "cmat2d.hpp"
#include "datetime/dtcalendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "geodesy/units.hpp"
#include "iers2010/doodson.hpp"
#include "stokes_coeffs.hpp"
#include <array>
#include <cstring>

namespace dso {
dso::iStatus parse_desai_ocean_pole_tide_coeffs(
    const char *fn, int max_degree, int max_order,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &Areal,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &Aimag,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &Breal,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &Bimag) noexcept;

/* @brief Fractional years since 01/01/2000
 * @param[in] mjdtt MJD of date of request [TT]
 * @return mjdtt - 01/01/2000 [TT] in fractional years
 */
inline double years_since_2000(const dso::TwoPartDate &mjdtt) noexcept {
  /* 01/01/2000 in MJD is: */
  const dso::TwoPartDate d2000(51544e0, 0e0);
  const auto dt_days = mjdtt - d2000;
  return dt_days.big() / 365.25e0 + dt_days.small() / 365.25e0;
}

/* @brief Secular pole coordinates (xp, yp) in [mas], computed as in IERS2010,
 *        Sec7.2.4 Eq. 21
 * @param[in] mjdtt MJD of date of request [TT]
 * @param[out] xs X component of secular pole coordinates, i.e. xp in [mas]
 * @param[out] ys Y component of secular pole coordinates, i.e. yp in [mas]
 */
inline int secular_pole(const dso::TwoPartDate &mjdtt, double &xs,
                        double &ys) noexcept {
  // secular pole (IERS2010, Sec7.2.4 Eq. 21)
  const double dt = years_since_2000(mjdtt);
  xs = 55e0 + 1.677e0 * dt;    // [mas]
  ys = 320.5e0 + 3.460e0 * dt; // [mas]
  return 0;
}

/* @brief Compute "wobble" variables m1 and m2 of polar motion (for pole
 *        tides)
 * m1 and m2 describe the time-dependent offset of the instantaneous rotation
 * pole from (xs, ys). See IERS2010, Sec7.2.4.
 * @param[in] xp Polar motion component X, in [as]
 * @param[in] yp Polar motion component Y, in [as]
 * @param[in] xs Secular pole component X, in [mas]
 * @param[in] ys Secular pole component Y, in [mas]
 * @param[out] m1 Woble X component in [as]
 * @param[out] m2 Woble Y component in [as]
 */
inline int wobble_components(double xp, double yp, double xs, double ys,
                             double &m1, double &m2) noexcept {
  /* "wobble" variables (IERS2010, Sec7.2.4 Eq. 25) */
  m1 = xp - xs * 1e-3;    // [as]
  m2 = -(yp - ys * 1e-3); // [as]
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
    int max_degree, int max_order, int num_header_lines,
    double scale = 1e0) noexcept;

class OceanTide {
private:
  std::vector<DoodsonOceanTideConstituent> doodsonFreqs;
  dso::StokesCoeffs dCS;
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> V; ///< workspace
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> W; ///< workspace
public:
  int max_degree() const noexcept { return dCS.max_degree(); }
  int max_order() const noexcept { return dCS.max_order(); }
  OceanTide(const std::vector<DoodsonOceanTideConstituent> &dfvec,
            int max_degree, int max_order, double GMearth = iers2010::GMe,
            double Rearth = iers2010::Re) noexcept
      : doodsonFreqs(dfvec), dCS(max_degree, max_order, GMearth, Rearth),
        V(max_degree + 3, max_degree + 3), W(max_degree + 3, max_degree + 3) {}

  int operator()(const dso::TwoPartDate &mjdtt, int max_degree,
                 int max_order) noexcept;
  int acceleration(const dso::TwoPartDate &mjdtt,
                   const Eigen::Matrix<double, 3, 1> &rsat,
                   Eigen::Matrix<double, 3, 1> &acc,
                   Eigen::Matrix<double, 3, 3> &acc_gradient,
                   int max_degree = -1, int max_order = -1) noexcept;
}; // OceanTide

class OceanPoleTide {
private:
  int max_degree{0};
  int max_order{0};
  /* Anm real part */
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> Ar;
  /* Anm imaginary part */
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> Ai;
  /* Bnm real part */
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> Br;
  /* Bnm imaginary part */
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> Bi;
  /* stokes coefficients (to be computed) */
  dso::StokesCoeffs dCS;
  /* degree-dependent coefficients for geopotential computation */
  // double *Rn{nullptr};
  std::vector<double> Rn;
  /* workspace */
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> V;
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> W;

public:
  OceanPoleTide(int maxdegree, int maxorder, const char *fn,
                double GMearth = iers2010::GMe, double Rearth = iers2010::Re);
  //~OceanPoleTide() noexcept;
  int operator()(const dso::TwoPartDate &mjdtt, double xp_sec, double yp_sec,
                 int maxdegre = -1, int maxorder = -1) noexcept;
  int acceleration(const dso::TwoPartDate &mjdtt, double xp_sec, double yp_sec,
                   const Eigen::Matrix<double, 3, 1> &rsat,
                   Eigen::Matrix<double, 3, 1> &acc,
                   Eigen::Matrix<double, 3, 3> &acc_gradient, int maxdegre = -1,
                   int maxorder = -1) noexcept;

}; /* OceanPoleTide */

class SolidEarthPoleTide {
  /* wobble variables in [asec] */
  double m1, m2;

  /* @brief Compute "wobble" variables m1 and m2 using Eqs. 24 and 25 from
   *        IERS2010, Sec7.2.4
   * The computed m1 and m2 values will be set in the instance's private
   * variables. Computed m1 and m2 are in [arcsec].
   *
   * @param[in] mjdtt Date of request in MJD, [TT]
   * @param[in] xp_sec Polar motion in [as] for the date of request
   * @param[in] yp_sec Polar motion in [as] for the date of request
   */
  int operator()(const dso::TwoPartDate &mjdtt, double xp_sec,
                 double yp_sec) noexcept {
    /* secular pole in [mas] */
    double xs, ys;
    secular_pole(mjdtt, xs, ys);
    /* wobble [as] */
    wobble_components(xp_sec, yp_sec, xs, ys, m1, m2);
    return 0;
  }

public:
  /* @brief Compute geopotential coefficient correction ΔC_21 and ΔS_21 due to
   *        Solid Earth pole tide, as descibed in IERS2010, Sec. 6.4
   * @param[in] mjdtt Date of request in MJD, [TT]
   * @param[in] xp_sec Polar motion in [as] for the date of request
   * @param[in] yp_sec Polar motion in [as] for the date of request
   * @return Effect of solid earth pole tide as correction in the normalized
   *         Stoke's coefficients C21 and S21. I.e. a pair is returned, as:
   *         [ΔC_21, ΔS_21]
   */
  auto delta_stokes_21(const dso::TwoPartDate &mjdtt, double xp_sec,
                       double yp_sec) noexcept {
    struct dCS21 {
      double dc21, ds21;
    };
    /* compute "wobble" variables of date */
    this->operator()(mjdtt, xp_sec, yp_sec);
    double dC21, dS21;
    /* compute solid earth pole tide */
    dC21 = -1.3331e-9 * (m1 + 0.0115e0 * m2);
    dS21 = -1.3331e-9 * (m2 - 0.0115e0 * m1);
    return dCS21{dC21, dS21};
  }

  int acceleration(const dso::TwoPartDate &mjdtt, double xp_sec, double yp_sec,
                   const Eigen::Matrix<double, 3, 1> &rsat,
                   Eigen::Matrix<double, 3, 1> &acc,
                   Eigen::Matrix<double, 3, 3> &acc_gradient,
                   double GMearth = iers2010::GMe,
                   double Rearth = iers2010::Re) noexcept;
}; /* SolidEarthPoleTide */

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
  int solid_earth_tide_step2_d(const dso::TwoPartDate &mjdtt, double &dC20,
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
  SolidEarthTide(double GMearth = iers2010::GMe, double Rearth = iers2010::Re,
                 double GMmoon = 0.49028010560e13,
                 double GMsun = 1.32712442076e20) noexcept
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

  /*int acceleration(const dso::TwoPartDate &mjdtt,
                   const dso::TwoPartDate &mjdut1,
                   double xp_sec, double yp_sec,
                   const Eigen::Matrix<double, 3, 1> &rsat,
                   const Eigen::Matrix<double, 3, 1> &rmoon,
                   const Eigen::Matrix<double, 3, 1> &rsun,
                   Eigen::Matrix<double, 3, 1> &acc,
                   Eigen::Matrix<double, 3, 3> &acc_gradient) noexcept;*/
}; // SolidEarthTide

} // namespace dso

#endif
