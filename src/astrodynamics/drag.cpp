#include "astrodynamics.hpp"
#include <iers2010/iersc.hpp>

/// I guess i don't really need the rbpn matrix ...
/// @brief Compute acceleration due to atmospheric drag
///
/// @param[in] rsat Satellite position in GCRS/ICRF [m]
/// @param[in] vsat Satellite velocity in GCRS/ICRF [m/sec]
/// @param[in] rbpn GCRs to true-of-date transformation matrix.
///            The matrix operates in the sense V(date) = rbpn * V(GCRS),
///            where the vector V(date) is with respect to the true equatorial
///            triad of date and the vector V(GCRS) is with respect to the
///            Geocentric Celestial Reference System (IAU, 2000).
///            See iers2010::iau::pnm06
/// @param[in] Area Cross-section [m^2]
/// @param[in] CD   Drag coefficient
/// @param[in] Mass Spacecraft mass [kg]
/// @param[in] Atmospheric density [kg/m^3]
Eigen::Matrix<double, 3, 1>
dso::drag_accel(const Eigen::Matrix<double, 3, 1> &rsat,
                const Eigen::Matrix<double, 3, 1> &vsat,
                const Eigen::Matrix<double, 3, 3> &rbpn, double Area, double CD,
                double Mass, double atmdens) noexcept {
  // earth angular velocity vector [rad/sec]
  constexpr const double omegav[] = {0e0, 0e0, iers2010::OmegaEarth};
  const Eigen::Matrix<double, 3, 1> omega{omegav};

  // position and velocity of satellite to true-of-date system
  const auto rtod = rbpn * rsat;
  const auto vtod = rbpn * vsat;

  // Velocity relative to the Earth's atmosphere
  const auto vrel = vtod - omega.cross(rtod);
  const auto vabs = vrel.norm();

  // Acceleration
  Eigen::Matrix<double, 3, 1> acc =
      -0.5e0 * CD * (Area / Mass) * atmdens * vabs * vrel;
  acc = rbpn.transpose() * acc;

  return acc;
}

Eigen::Matrix<double, 3, 1>
dso::drag_accel(const Eigen::Matrix<double, 3, 1> &rsat,
                const Eigen::Matrix<double, 3, 1> &vsat, double Area, double CD,
                double Mass, double atmdens) noexcept {
  // earth angular velocity vector [rad/sec]
  constexpr const double omegav[] = {0e0, 0e0, iers2010::OmegaEarth};
  const Eigen::Matrix<double, 3, 1> omega{omegav};

  // Velocity relative to the Earth's atmosphere
  const Eigen::Matrix<double, 3, 1> vrel = vsat - omega.cross(rsat);
  // const auto vabs = vrel.norm();

  // Acceleration
  Eigen::Matrix<double, 3, 1> acc = -0.5e0 * CD * (Area / Mass) * atmdens *
                                    vrel.squaredNorm() * vrel.normalized();

  return acc;
}

Eigen::Matrix<double, 3, 1> dso::drag_accel(
    const Eigen::Matrix<double, 3, 1> &rsat,
    const Eigen::Matrix<double, 3, 1> &vsat, double Area, double CD,
    double Mass, double atmdens, const Eigen::Matrix<double, 3, 1> &datmdensdr,
    Eigen::Matrix<double, 3, 3> &daccdr, Eigen::Matrix<double, 3, 3> &daccdv,
    Eigen::Matrix<double, 3, 1> &daccdC) noexcept {
  // earth angular velocity vector [rad/sec]
  constexpr const double omegav[] = {0e0, 0e0, iers2010::OmegaEarth};
  const Eigen::Matrix<double, 3, 1> omega{omegav};

  // Velocity relative to the Earth's atmosphere
  const Eigen::Matrix<double, 3, 1> vrel = vsat - omega.cross(rsat);
  // const auto vabs = vrel.norm();

  // Acceleration
  Eigen::Matrix<double, 3, 1> acc = -0.5e0 * CD * (Area / Mass) * atmdens *
                                    vrel.squaredNorm() * vrel.normalized();

  // partials w.r.t drag coefficient (Montenbruck, 7.80)
  daccdC = -0.5e0 * (Area / Mass) * atmdens * vrel.norm() * vrel;

  // partials w.r.t velocity (Montenbruck, 7.81)
  daccdv = -0.5e0 * CD * (Area / Mass) * atmdens *
           (vrel.norm() * Eigen::Matrix<double, 3, 3>::Identity() +
            vrel * vrel.transpose() / vrel.norm());

  // Partials w.r.t. position vector (Montenbruck, 7.84)
  const double _data[] = {
      0e0, iers2010::OmegaEarth, 0e0, -iers2010::OmegaEarth, 0e0, 0e0, 0e0, 0e0,
      0e0};
  const Eigen::Matrix<double, 3, 3> XOmega(_data);
  daccdr = -0.5e0 * CD * (Area / Mass) * vrel.norm() * vrel *
               datmdensdr.transpose() -
           daccdv * XOmega;

  return acc;
}

Eigen::Matrix<double, 3, 1>
dso::drag_accel(const Eigen::Matrix<double, 3, 1> &vrel, double AM, double Cd,
                double atmdens, const Eigen::Matrix<double, 3, 1> &datmdensdr,
                Eigen::Matrix<double, 3, 3> &daccdr,
                Eigen::Matrix<double, 3, 3> &daccdv,
                Eigen::Matrix<double, 3, 1> &daccdC) noexcept {

  // vrel unit vector
  // const Eigen::Matrix<double, 3, 1> erel = vrel.normalized();
  // vrel norm
  const double Vrel = vrel.norm();
  
  // Acceleration
  Eigen::Matrix<double, 3, 1> acc = -0.5e0 * Cd * AM * atmdens * Vrel * vrel;

  // partials w.r.t drag coefficient (Montenbruck, 7.80)
  daccdC = -0.5e0 * AM * atmdens * Vrel * vrel;

  // partials w.r.t velocity (Montenbruck, 7.81)
  daccdv = -0.5e0 * Cd * AM * atmdens *
           (Vrel * Eigen::Matrix<double, 3, 3>::Identity() +
            vrel * vrel.transpose() / Vrel);

  // Partials w.r.t. position vector (Montenbruck, 7.84)
  const double _data[] = {
      0e0, iers2010::OmegaEarth, 0e0, -iers2010::OmegaEarth, 0e0, 0e0, 0e0, 0e0,
      0e0};
  const Eigen::Matrix<double, 3, 3> XOmega(_data);
  daccdr =
      -0.5e0 * Cd * AM * Vrel * vrel * datmdensdr.transpose() - daccdv * XOmega;

  return acc;
}
