#ifndef __ORBIT_INTEGRATION_PARAMETERS_HPP__
#define __ORBIT_INTEGRATION_PARAMETERS_HPP__

#include "atmosphere.hpp"
#include "atmosphere/dtm2020/dtm2020.hpp"
#include "egravity.hpp"
#include "eigen3/Eigen/Eigen"
#include "iers2010/eop.hpp"
#include "tides.hpp"
#include "planetpos.hpp"
#include "satellites.hpp"
#include "satellites/jason3_quaternions.hpp"
#include "attitude.hpp"
#include <cassert>

namespace dso {

/* @brief A structure to hold orbit integration parameters; it is a
 *        collection of data and parameters to be used in the computation
 *        of (the system of) variational equations.
 */
class IntegrationParameters {
  /* time in TAI */
  dso::TwoPartDate mjd_tai;
  
  /* EOP parameters Look-up table */
  const dso::EopLookUpTable &eopLUT;

  /* Earth gravity field */
  dso::EarthGravity *egravity;
  
  /* Sun/Moon gravitational parameters, in [m^3/ sec^2] */
  double GMSun, GMMon;

  /* Ocean Tides */
  dso::OceanTide *octide{nullptr};

  /* Earth Tides */
  dso::SolidEarthTide *setide{nullptr};
  
  /* Pole tides */
  dso::SolidEarthPoleTide *psetide{nullptr};
  dso::OceanPoleTide *poctide{nullptr};
  
  /* Satellite-specific information */
  SvFrame *svFrame{nullptr};

  /* atmospheric density model and data feed */
  dso::Dtm2020 Dtm20;

  /* setup dynamic parameters */
  double Cd = 2e0;
  double Cr = 1.5e0;

  /* workspace -- these should match the maximum degree and order */
  dso::StokesCoeffs geopotential;
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> V,W;

public:
  IntegrationParameters(const dso::EopLookUpTable &eoptable_,
                        dso::EarthGravity *earth_gravity,
                        dso::OceanTide *ocean_tide,
                        dso::SolidEarthTide *solid_earth_tide,
                        dso::SolidEarthPoleTide *solid_earth_pole_tide,
                        dso::OceanPoleTide *ocean_pole_tide,
                        const char *pck_kernel,
                        const char *dtm2020datafile);

  dso::TwoPartDate &reference_epoch() noexcept { return mjd_tai; }
  dso::TwoPartDate reference_epoch() const  noexcept { return mjd_tai; }
  const dso::EopLookUpTable &eop_lookup_table() const noexcept { return eopLUT; }
  double &drag_ceofficient() noexcept {return Cd;}
  double &srp_ceofficient() noexcept {return Cr;}
  dso::EarthGravity *earth_gravity() noexcept {return egravity;};
  dso::OceanTide *ocean_tide() noexcept {return octide;}
  dso::SolidEarthTide *solid_earth_tide() noexcept {return setide;}
  dso::SolidEarthPoleTide *solid_earth_pole_tide() noexcept {return psetide;}
  dso::OceanPoleTide *ocean_pole_tide() noexcept {return poctide;}
  dso::TwoPartDate tai() const noexcept {return mjd_tai;}
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &workspace_v() noexcept {return V;}
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &workspace_w() noexcept {return W;}

  dso::StokesCoeffs &accumulate_geopotential_coeffs() noexcept {
    geopotential.clear();
    if (egravity) geopotential+=egravity->geopotential_coeffs();
    if (octide) geopotential+=octide->geopotential_coeffs();
    if (setide) geopotential+=setide->geopotential_coeffs();
    if (poctide) geopotential+=poctide->geopotential_coeffs();
    if (psetide) {
      geopotential.C(2,1) += psetide->dC21();
      geopotential.S(2,1) += psetide->dS21();
    }
    return geopotential;
  }

  void set_flux_data(dso::SolarActivityData &data) noexcept {
    Dtm20.set_flux_data(data);
  }

  SvFrame* svframe() noexcept {return svFrame;}

  void set_sv_frame(Eigen::Matrix<double, 3, 1> vcog,
                    Eigen::Matrix<double, 3, 1> varp, const char *qfn, const char *sfn,
                    double mass) {
    svFrame = new dso::SvFrame(vcog,varp,qfn,sfn,mass);
  }

}; /* Integration Parameters */

/// @brief Comnpute third-body, Sun- and Moon- induced acceleration on an
///        orbiting satellite.
/// @warning Note that the function asserts that the respectice SPICE kernels
///        are already loaded (it will call dso::cspice::j2planet_pos_from
///        function).
/// @note The gravitational constants of Sun and Moon, should be retrieved by
///        a corresponding SPICE kernel. See dso::getSunMoonGM.
/// @param[in] mjd_tai TAI date as MJD
/// @param[in] GMSun Gravitational parameter for Sun (see Note 1)
/// @param[in] GMMon Gravitational parameter for Moon (see Note 1)
/// @param[in] rsat Position vector of the satellite in GCRS [m]
/// @param[out] sun_acc Sun-induced acceleration on satellite [m/sec^2]
/// @param[out] mon_acc Moon-induced acceleration on satellite [m/sec^2]
/// @param[out] sun_pos Sun position in J2000 [m]
/// @param[out] mon_partials Partials of the Moon-induced acceleration w.r.t
///       satellite position vecot (=r), aka d(acc)/dr
void SunMoon(const dso::TwoPartDate &mjd_tai, const Eigen::Matrix<double, 3, 1> &rsat,
             double GMSun, double GMMon, Eigen::Matrix<double, 3, 1> &sun_acc,
             Eigen::Matrix<double, 3, 1> &mon_acc,
             Eigen::Matrix<double, 3, 1> &sun_pos,
             Eigen::Matrix<double, 3, 1> &mon_pos,
             Eigen::Matrix<double, 3, 3> &gradient) noexcept;

void VariationalEquations(double tsec_away, const Eigen::VectorXd &yPhi,
                          Eigen::Ref<Eigen::VectorXd> yPhiP,
                          dso::IntegrationParameters *params) noexcept;

/*
void VariationalEquations(double tsec_away, const Eigen::VectorXd &yPhi,
                          Eigen::Ref<Eigen::VectorXd> yPhiP,
                          dso::IntegrationParameters &params) noexcept;
void VariationalEquations_mg(double tsec_away, const Eigen::VectorXd &yPhi,
                          Eigen::Ref<Eigen::VectorXd> yPhiP,
                          dso::IntegrationParameters &params) noexcept;

void VariationalEquations_ta(double tsec_away, const Eigen::VectorXd &yPhi,
                          Eigen::Ref<Eigen::VectorXd> yPhiP,
                          dso::IntegrationParameters &params) noexcept;
void VariationalEquations_thread(double tsec_away, const Eigen::VectorXd &yPhi,
                          Eigen::Ref<Eigen::VectorXd> yPhiP,
                          dso::IntegrationParameters &params) noexcept;
void noVariationalEquations(double tsec_away, const Eigen::VectorXd &yPhi,
                          Eigen::Ref<Eigen::VectorXd> yPhiP,
                          dso::IntegrationParameters &params) noexcept;
*/
} // namespace dso

#endif
