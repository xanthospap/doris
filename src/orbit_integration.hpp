#ifndef __ORBIT_INTEGRATION_PARAMETERS_HPP__
#define __ORBIT_INTEGRATION_PARAMETERS_HPP__

#include "egravity.hpp"
#include "eop.hpp"
#include "planetpos.hpp"
#include "satellites.hpp"
#include "satellites/jason3_quaternions.hpp"
#include "atmosphere.hpp"
#include "eigen3/Eigen/Eigen"
#include <cassert>

namespace dso {

/// @brief A structure to hold orbit integration parameters; it is a
///        collection of data and parameters to be used in the computation
///        of (the system of) variational equations.
struct IntegrationParameters {
  ///< time in TAI
  double mjd_tai;
  ///< EOP parameters Look-up table
  const dso::EopLookUpTable &eopLUT;
  ///< gravity harmonics
  const dso::HarmonicCoeffs &harmonics;
  ///< memmory for Lagrange polynomials
  dso::Mat2D<dso::MatrixStorageType::Trapezoid> *Lagrange_V{nullptr};
  dso::Mat2D<dso::MatrixStorageType::Trapezoid> *Lagrange_W{nullptr};
  ///< degree and order of geopotential harmonics
  int degree, order;
  ///< Sun/Moon gravitational parameters
  double GMSun, GMMon;
  ///< Satellite Macromodel and number of individual flat plates
  const MacroModelComponent *macromodel{nullptr};
  int numMacroModelComponents{0};
  dso::JasonQuaternionHunter *qhunt{nullptr};
  const double *SatMass{nullptr};
  /// Drag-related stuff
  dso::nrlmsise00::InParams<
      dso::nrlmsise00::detail::FluxDataFeedType::ST_CSV_SW> *AtmDataFeed;
  dso::Nrlmsise00 *nrlmsise00;
  //Eigen::Matrix<double,3,1> ddragdC;
  Eigen::VectorXd *estimates;

  IntegrationParameters(int degree_, int order_,
                        const dso::EopLookUpTable &eoptable_,
                        const dso::HarmonicCoeffs &harmonics_,
                        const char *pck_kernel) noexcept
      : eopLUT(eoptable_), harmonics(harmonics_),
        Lagrange_V{new dso::Mat2D<dso::MatrixStorageType::Trapezoid>(
            degree_ + 3, order_ + 3)},
        Lagrange_W{new dso::Mat2D<dso::MatrixStorageType::Trapezoid>(
            degree_ + 3, order_ + 3)},
        degree(degree_), order(order_) {
    assert(degree_ == harmonics_.degree());
    // gravitational parameters
    assert(!dso::get_sun_moon_GM(pck_kernel, GMSun, GMMon));
  };

  ~IntegrationParameters() noexcept {
    if (Lagrange_V)
      delete Lagrange_V;
    if (Lagrange_W)
      delete Lagrange_W;
  }
}; // Integration Parameters

/// @brief Compute the Terrestrial-to-Celestial (aka ITRS to GCRS) matrix
///        using IAU2000 model
/// @param[in] mjd_tai TAI date as MJD
/// @param[in] eop_table A dso::EopLookUpTable to be used to get EOP
///            parameters for the given date
/// @param[out] ditrs2gcrs Derivative of the returned Terrestrial-to-Celestial
///            matrix w.r.t time. For the computation of the derivative,
///            assuming (Petit et al, 2010): [GCRS] = Q(t) R(t) W(t) [ITRS]
///            we consider the Q and W matrices as constant, and only compute
///            the derivative of the R matrix. That is (Montenbruck et al, 2012,
///            Chapter 5, Excercise 5.2):
///            if U = Q(t) * R(t) * W(t),
///            dU/dt = Q(t) * dR(t)/dt * W(t) and dR(t)/dt = ω * T * R(t)
///                     | 0 +1 0 |
///            with T = |-1  0 0 |
///                     | 0  0 0 |
///            (note that the above derivation follows the Celestial-to-
///            Terrestrial matrix; we want its transpose).
Eigen::Matrix<double, 3, 3>
itrs2gcrs(double mjd_tai, const dso::EopLookUpTable &eop_table,
          Eigen::Matrix<double, 3, 3> &ditrs2gcrs) noexcept;

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
void SunMoon(double mjd_tai, const Eigen::Matrix<double, 3, 1> &rsat,
             double GMSun, double GMMon, Eigen::Matrix<double, 3, 1> &sun_acc,
             Eigen::Matrix<double, 3, 1> &mon_acc,
             Eigen::Matrix<double, 3, 1> &sun_pos,
             Eigen::Matrix<double, 3, 3> &mon_partials) noexcept;

void VariationalEquations(double tsec, const Eigen::VectorXd &yPhi,
                          Eigen::Ref<Eigen::VectorXd> yPhiP,
                          dso::IntegrationParameters &params) noexcept;

} // namespace dso

#endif
