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
#include <iers2010/iersc.hpp>

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

/* @brief Ocean Tide information for an individual wave
 * This class holds information, i.e. DelC+, DelS+, DelC-, DelS- coefficients
 * up to degree 'maxl' and order 'maxm', for a given tidal wave, identified by
 * its Doodson number 'doodson' (see Eq. 6.15 in IERS2010 standards).
 * Coefficients are usually read off from a relevant input file.
 */
class DoodsonOceanTideConstituent {
private:
  /* Doodson number of the wave */
  DoodsonNumber doodson;
  /* max degree of coefficients */
  int maxl;
  /* max order of coefficients */
  int maxm;
  /* coefficients DelC+, i.e. C_nm^(+) */
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> DelCpl;
  /* coefficients DelS+, i.e. S_nm^(+) */
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> DelSpl;
  /* coefficients DelC-, i.e. C_nm^(-) */
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> DelCmi;
  /* coefficients DelS-, i.e. S_nm^(-) */
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> DelSmi;

public:
  DoodsonOceanTideConstituent(const dso::DoodsonNumber d, int max_degree,
                              int max_order);

  int max_degree() const noexcept { return maxl; }
  int max_order() const noexcept { return maxm; }
  const DoodsonNumber &doodson_number() const noexcept { return doodson; }
  const double &delCp(int l, int m) const noexcept {
    return DelCpl.operator()(l, m);
  }
  const double &delSp(int l, int m) const noexcept {
    return DelSpl.operator()(l, m);
  }
  const double &delCm(int l, int m) const noexcept {
    return DelCmi.operator()(l, m);
  }
  const double &delSm(int l, int m) const noexcept {
    return DelSmi.operator()(l, m);
  }
  double &delCp(int l, int m) noexcept { return DelCpl.operator()(l, m); }
  double &delSp(int l, int m) noexcept { return DelSpl.operator()(l, m); }
  double &delCm(int l, int m) noexcept { return DelCmi.operator()(l, m); }
  double &delSm(int l, int m) noexcept { return DelSmi.operator()(l, m); }

  /* resize instance to hold coefficients up to maxDegree (inclusive) */
  void resize(int maxDegree) noexcept;

  /* set all of C_nm^(+) S_nm^(+) C_nm^(-) and S_nm^(-) coeffs to zero */
  void clear_coefficients() noexcept;
}; /* DoodsonOceanTideConstituent */

iStatus memmap_octide_coefficients(
    const char *fn, std::vector<dso::DoodsonOceanTideConstituent> &freqs,
    int max_degree, int max_order, double scale = 1e0) noexcept;

/* @brief Match given waves from an DoodsonOceanTideConstituent vector
 * @param[in] waves A vector of DoodsonNumber of interest
 * @param[in] freqs Original vector of DoodsonOceanTideConstituent (to be
 *            filtered)
 * @return std::vector<dso::DoodsonOceanTideConstituent> A vector of
 *         DoodsonOceanTideConstituent, containing only waves included both
 *         in waves and in freqs vectors.
 */
std::vector<dso::DoodsonOceanTideConstituent> filter_waves(
    const std::vector<DoodsonNumber> &waves,
    const std::vector<dso::DoodsonOceanTideConstituent> &freqs) noexcept;

/* @brief A class to hold and compute ocean tidal disturbances
 * This class i based on model described in IERS-2010 standards, Sec. 6.3.
 * It holds a list of 'main' tidal waves, along with their coefficients i.e.
 * C_nm^(+) S_nm^(+) C_nm^(-) and S_nm^(-) up to given degree and order (n,m).
 * For each such wave, the relevant information are kept in an individual
 * DoodsonOceanTideConstituent instance.
 * Computation of geopotential ΔCnm and ΔSnm components is based on Eq. 6.15
 */
class OceanTide {
private:
  /* hold one DoodsonOceanTideConstituent for each of the main waves of the
   * model. This instance is usually read-off of from an input file (e.g.
   * ftp://tai.bipm.org/iers/conv2010/chapter6/tidemodels/fes2004_Cnm-Snm.dat
   * and function memmap_octide_coefficients)
   */
  std::vector<DoodsonOceanTideConstituent> doodsonFreqs;
  /* ΔCnm and ΔSnm after computation */
  dso::StokesCoeffs dCS;

public:
  int max_degree() const noexcept { return dCS.max_degree(); }
  int max_order() const noexcept { return dCS.max_order(); }
  int num_main_waves() const noexcept { return doodsonFreqs.size(); }
  dso::StokesCoeffs &geopotential_coeffs() noexcept {return dCS;}

  OceanTide(const char *fn, double scale, int max_degree, int max_order,
            double GMearth = iers2010::GMe, double Rearth = iers2010::Re);

  /* Compute geopotential coefficients ΔC and ΔS and store them in dCS */
  iStatus operator()(const dso::TwoPartDate &mjdtt,
                     const dso::TwoPartDate &mjdut1, int max_degree = -1,
                     int max_order = -1) noexcept;
  
  /* Compute geopotential acceleration at given point. 
   *
   * If you use this function a lot, then you should pre-allocate the V and W
   * matrices (with size >= degree + 2) and pass in the respective pointers.
   * Otherwise, the function will try to allocate storage each time it is
   * called.
   * The function will call the operator() on the instance to compute the 
   * geopotential coefficients, and then compute acceleration (and gradient) 
   * based on an SH expansion.
   *
   * All units is SI. Reference frame is ITRF.
   */
  iStatus acceleration(
      const dso::TwoPartDate &mjdtt, const dso::TwoPartDate &mjdut1,
      const Eigen::Matrix<double, 3, 1> &rsat, Eigen::Matrix<double, 3, 1> &acc,
      Eigen::Matrix<double, 3, 3> &acc_gradient, int max_degree = -1,
      int max_order = -1,
      dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *V = nullptr,
      dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *W =
          nullptr) noexcept;
}; /* OceanTide */

class AtmosphericTide {
private:
  /* hold one DoodsonOceanTideConstituent for each of the main waves of the
   * model. This instance is usually read-off of from an input file (e.g.
   * atmosTides_AOD1BRL06.potential.iers.txt
   * and function memmap_octide_coefficients)
   */
  std::vector<DoodsonOceanTideConstituent> doodsonFreqs;
  /* ΔCnm and ΔSnm after computation */
  dso::StokesCoeffs dCS;
  /* workspace matrix */
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> V;
  /* workspace matrix */
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> W;

public:
  int max_degree() const noexcept { return dCS.max_degree(); }
  int max_order() const noexcept { return dCS.max_order(); }
  int num_main_waves() const noexcept { return doodsonFreqs.size(); }

  void filter_waves(const std::vector<DoodsonNumber> &waves) noexcept {
    doodsonFreqs = dso::filter_waves(waves, doodsonFreqs);
  }

  AtmosphericTide(const char *fn, double scale, int max_degree, int max_order,
                  double GMearth = iers2010::GMe, double Rearth = iers2010::Re);

  iStatus operator()(const dso::TwoPartDate &mjdtt,
                     const dso::TwoPartDate &mjdut1, int max_degree = -1,
                     int max_order = -1) noexcept;
  iStatus acceleration(const dso::TwoPartDate &mjdtt,
                       const dso::TwoPartDate &mjdut1,
                       const Eigen::Matrix<double, 3, 1> &rsat,
                       Eigen::Matrix<double, 3, 1> &acc,
                       Eigen::Matrix<double, 3, 3> &acc_gradient,
                       int max_degree = -1, int max_order = -1) noexcept;
}; /* AtmosphericTide */

class OceanPoleTide {
private:
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
  std::vector<double> Rn;

public:
  int max_degree() const noexcept { return dCS.max_degree(); }
  int max_order() const noexcept { return dCS.max_order(); }

  OceanPoleTide(int maxdegree, int maxorder, const char *fn,
                double GMearth = iers2010::GMe, double Rearth = iers2010::Re);

  dso::StokesCoeffs &geopotential_coeffs() noexcept {return dCS;}

  /* Compute geopotential coefficients ΔC and ΔS and store them in dCS */
  iStatus operator()(const dso::TwoPartDate &mjdtt, double xp_sec,
                     double yp_sec, int maxdegre = -1,
                     int maxorder = -1) noexcept;

  /* Compute geopotential acceleration at given point.
   *
   * If you use this function a lot, then you should pre-allocate the V and W
   * matrices (with size >= degree + 2) and pass in the respective pointers.
   * Otherwise, the function will try to allocate storage each time it is
   * called.
   * The function will call the operator() on the instance to compute the
   * geopotential coefficients, and then compute acceleration (and gradient)
   * based on an SH expansion.
   *
   * All units is SI. Reference frame is ITRF.
   */
  iStatus acceleration(
      const dso::TwoPartDate &mjdtt, double xp_sec, double yp_sec,
      const Eigen::Matrix<double, 3, 1> &rsat, Eigen::Matrix<double, 3, 1> &acc,
      Eigen::Matrix<double, 3, 3> &acc_gradient, int maxdegre = -1,
      int maxorder = -1,
      dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *V = nullptr,
      dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *W =
          nullptr) noexcept;

}; /* OceanPoleTide */

class SolidEarthPoleTide {
  /* wobble variables in [asec] */
  double m1, m2;
  /* geopotential coefficients (after computation) */
  double mdC21, mdS21;

  /* @brief Compute "wobble" variables m1 and m2 using Eqs. 24 and 25 from
   *        IERS2010, Sec7.2.4
   * The computed m1 and m2 values will be set in the instance's private
   * variables. Computed m1 and m2 are in [arcsec].
   *
   * @param[in] mjdtt Date of request in MJD, [TT]
   * @param[in] xp_sec Polar motion in [as] for the date of request
   * @param[in] yp_sec Polar motion in [as] for the date of request
   */
  int wobble(const dso::TwoPartDate &mjdtt, double xp_sec,
                 double yp_sec) noexcept {
    /* secular pole in [mas] */
    double xs, ys;
    secular_pole(mjdtt, xs, ys);
    /* wobble [as] */
    wobble_components(xp_sec, yp_sec, xs, ys, m1, m2);
    return 0;
  }

public:
  constexpr int max_degree() const noexcept {return 2;}
  constexpr int max_order() const noexcept {return 1;}
  double dC21() const noexcept {return mdC21;}
  double dS21() const noexcept {return mdS21;}

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
    this->wobble(mjdtt, xp_sec, yp_sec);
    /* compute solid earth pole tide */
    mdC21 = -1.3331e-9 * (m1 + 0.0115e0 * m2);
    mdS21 = -1.3331e-9 * (m2 - 0.0115e0 * m1);
    return dCS21{mdC21, mdS21};
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
  /* gravitational constants for Sun and Moon in SI */
  const double GM_moon, GM_sun;
  /* geopotential, i.e. Stoke's coefficients (n,m) = (4,4) */
  dso::StokesCoeffs cs;

  int solid_earth_tide_step1(const Eigen::Matrix<double, 3, 1> &rmoon,
                             const Eigen::Matrix<double, 3, 1> &rsun,
                             std::array<double, 12> &dC,
                             std::array<double, 12> &dS) noexcept;

  int solid_earth_tide_step2(const dso::TwoPartDate &mjdtt,
                             const dso::TwoPartDate &mjdut1, double &dC20,
                             double &dC21, double &dS21, double &dC22,
                             double &dS22) const noexcept;

public:
  int max_degree() const noexcept {return cs.max_degree(); }
  int max_order() const noexcept {return cs.max_order(); }

  /* @brief Constructor
   * @param GMearth Standard gravitational parameter μ=GM for the Earth
   * [m^2/s^2]
   * @param Rearth Equatorial radius of Earth [m]
   * @param GMmoon Standard gravitational parameter μ=GM for the Moon
   *        [m^2/s^2]
   * @param GMsun Standard gravitational parameter μ=GM for the Moon
   *        [m^2/s^2]
   */
  SolidEarthTide(double GMearth = iers2010::GMe, double Rearth = iers2010::Re,
                 double GMmoon = 0.49028010560e13,
                 double GMsun = iers2010::GMSun) noexcept
      : GM_moon(GMmoon), GM_sun(GMsun), cs(degree, degree, GMearth, Rearth)
  {}

  dso::StokesCoeffs &geopotential_coeffs() noexcept {return cs;}

  /* Compute geopotential coefficients ΔC and ΔS and store them in cs */
  int operator()(const dso::TwoPartDate &mjdtt, const dso::TwoPartDate &mjdut1,
                 const Eigen::Matrix<double, 3, 1> &rmoon,
                 const Eigen::Matrix<double, 3, 1> &rsun) noexcept;

  /* Compute geopotential acceleration at given point. 
   *
   * If you use this function a lot, then you should pre-allocate the V and W
   * matrices (with size >= degree + 2) and pass in the respective pointers.
   * Otherwise, the function will try to allocate storage each time it is
   * called.
   * The function will call the operator() on the instance to compute the 
   * geopotential coefficients, and then compute acceleration (and gradient) 
   * based on an SH expansion.
   *
   * All units is SI. Reference frame is ITRF.
   */
  dso::iStatus acceleration(
      const dso::TwoPartDate &mjdtt, const dso::TwoPartDate &mjdut1,
      const Eigen::Matrix<double, 3, 1> &rsat,
      const Eigen::Matrix<double, 3, 1> &rmoon,
      const Eigen::Matrix<double, 3, 1> &rsun, Eigen::Matrix<double, 3, 1> &acc,
      Eigen::Matrix<double, 3, 3> &acc_gradient,
      dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *V = nullptr,
      dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *W =
          nullptr) noexcept;
}; /* SolidEarthTide */

} /* namespace dso */

#endif
