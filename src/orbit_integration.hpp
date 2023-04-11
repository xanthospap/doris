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
#include <cassert>
#include <datetime/dtcalendar.hpp>

namespace dso {

class SvFrame {
  /* CoG point in SV-fixed RF */
  Eigen::Matrix<double, 3, 1> cog_sf;
  /* Antenna RP point in SV-fixed RF */
  Eigen::Matrix<double, 3, 1> arp_sf;
  /* quaternion hunter for attitude */
  // dso::JasonQuaternionHunter qhunt;
  dso::JasonAttitude jatt;
  /* satellite macromodel */
  dso::MacroModel<SATELLITE::Jason3> mplates;
  double sat_mass{0e0};

public:
  SvFrame(Eigen::Matrix<double, 3, 1> vcog, Eigen::Matrix<double, 3, 1> varp,
          const char *qfn, const char *bfn, double mass)
      : cog_sf(vcog), arp_sf(varp), jatt(qfn,bfn), sat_mass(mass) {}

  double mass() const noexcept {return sat_mass;}
  const dso::MacroModel<SATELLITE::Jason3> &plates() const noexcept {return mplates;}

  int arp2cog(const dso::TwoPartDate &tai, Eigen::Matrix<double, 3, 1> &dr) {
    /* r_sat_fixed = Q * r_sat_gcrs */
    Eigen::Quaternion<double> q;
    if (jatt.get_at(tai, q)) {
      fprintf(stderr, "[ERROR] Failed to find quaternion for datetime\n");
      return 1;
    }
    dr = q.conjugate().normalized() * (cog_sf - arp_sf);
    return 0;
  }

  /*
  int get_attitude_quaternion(const dso::TwoPartDate &tai, Eigen::Quaternion<double> &q) noexcept {
    if (qhunt.get_at(tai, q)) {
      fprintf(stderr, "[ERROR] Failed to find quaternion for datetime\n");
      return 1;
    }
    return 0;
  }
  */

  /* returns (S / m) where S is the projected area */
  int projected_area(const dso::TwoPartDate &tai,
                     const Eigen::Matrix<double, 3, 1> &vrel, double &b) noexcept {
    double left, right;
    /* attitude matrix */
    Eigen::Quaternion<double> q;
    if (jatt.get_at(tai, q, left, right)) {
      fprintf(stderr, "[ERROR] Failed to find quaternion for datetime\n");
      return 1;
    }
    
    /* normalized relative velocity */
    const auto v = vrel.normalized();
    
    /* iterate through SV plates, no solar array yet */
    double S = 0e0;
    const int numPlates = mplates.NumPlates;
    const MacroModelComponent *plate = nullptr;
    for (int i=0; i<numPlates-2; i++) {
      plate = mplates.mmcomponents + i;
      /* (unit) vector of plate, normal in sv-fixed RF */
      Eigen::Matrix<double,3,1> n(plate->m_normal);
      const double cosA = (q.conjugate().normalized()*n).transpose() * v;
      if (cosA > 0e0) {
        S += (plate->m_surf * cosA);
      }
    }

    int i;
    /* left panel */
    {
    i = numPlates - 2;
    plate = mplates.mmcomponents + i;
    /* (unit) vector of plate, normal in sv-fixed RF */
    const Eigen::Matrix<double,3,3> R = Eigen::AngleAxis<double>(left, -Eigen::Vector3d::UnitY()).toRotationMatrix().transpose();
    Eigen::Matrix<double, 3, 1> n(R * Eigen::Matrix<double,3,1>{plate->m_normal});
    const double cosA = (q.conjugate().normalized() * n).transpose() * v;
    if (cosA > 0e0) {
      S += (plate->m_surf * cosA);
    }
    }

    {
    ++i;
    plate = mplates.mmcomponents + i;
    /* (unit) vector of plate, normal in sv-fixed RF */
    const auto R = Eigen::AngleAxis<double>(right, Eigen::Matrix<double, 3, 1>{0e0, 1e0, 0e0});
    Eigen::Matrix<double, 3, 1> n(
            R.toRotationMatrix()
            .transpose() *
        Eigen::Matrix<double,3,1>{plate->m_normal});
    const double cosA = (q.conjugate().normalized() * n).transpose() * v;
    if (cosA > 0e0) {
      S += (plate->m_surf * cosA);
    }
    }

    /* compute ballistic coefficient, without the Cd */
    b = S/sat_mass;
    return 0;
  }

};


/// @brief A structure to hold orbit integration parameters; it is a
///        collection of data and parameters to be used in the computation
///        of (the system of) variational equations.
struct IntegrationParameters {
  ///< time in TAI
  dso::TwoPartDate mjd_tai;
  dso::TwoPartDate &reference_epoch() noexcept { return mjd_tai; }
  dso::TwoPartDate reference_epoch() const  noexcept { return mjd_tai; }
  ///< EOP parameters Look-up table
  const dso::EopLookUpTable &eopLUT;
  const dso::EopLookUpTable &eop_lookup_table() const noexcept { return eopLUT; }
  ///< gravity harmonics
  const dso::HarmonicCoeffs &harmonics;
  ///< degree and order of geopotential harmonics
  int degree, order;
  ///< Sun/Moon gravitational parameters, in [km^3/ sec^2]
  double GMSun, GMMon;
  ///< Ocean Tides
  dso::OceanTide *octide;
  ///< Earth Tides
  dso::SolidEarthTide *setide;
  ///< Satellite-specific information
  SvFrame *svFrame{nullptr};
  /* atmospheric density model and data feed */
  dso::Dtm2020 Dtm20;
  /* setup dynamic parameters */
  double Cd = 2e0;
  double &drag_ceofficient() noexcept {return Cd;}
  double Cr = 1.5e0;
  double &srp_ceofficient() noexcept {return Cr;}
  ///< memmory for Lagrange polynomials
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *V{nullptr};
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *W{nullptr};
  
  IntegrationParameters(int degree_, int order_,
                        const dso::EopLookUpTable &eoptable_,
                        const dso::HarmonicCoeffs &harmonics_,
                        const char *pck_kernel, const char *dtm2020datafile) noexcept
      : eopLUT(eoptable_), harmonics(harmonics_), degree(degree_),
        order(order_), Dtm20(dtm2020datafile),
        V(new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
            degree_ + 3, degree_ + 3)),
        W(new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
            degree_ + 3, degree_ + 3)) {
    assert(degree_ == harmonics_.max_degree());
    assert(!dso::get_sun_moon_GM(pck_kernel, GMSun, GMMon)); // [km^3/ sec^2]
  };

  void set_sv_frame(Eigen::Matrix<double, 3, 1> vcog,
                    Eigen::Matrix<double, 3, 1> varp, const char *qfn,
                    const char *bfn, double mass) {
    svFrame = new dso::SvFrame(vcog, varp, qfn, bfn, mass);
  }

  // TODO
  ~IntegrationParameters() noexcept {
    if (V)
      delete V;
    if (W)
      delete W;
  }
}; // Integration Parameters

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
                          dso::IntegrationParameters &params) noexcept;
void VariationalEquations2(double tsec_away, const Eigen::VectorXd &yPhi,
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
} // namespace dso

#endif
