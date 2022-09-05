#ifndef __IDS_CRYOSAT2_SATELLITE_DETAILS_HPP__
#define __IDS_CRYOSAT2_SATELLITE_DETAILS_HPP__

#include "datetime/dtcalendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "satellites.hpp"
#include <cassert>

namespace dso {
template <> struct SatelliteInfo<SATELLITE::Cryosat2> {
  /// @ref https://cddis.nasa.gov/Techniques/sp3c_satlist.html
  static constexpr char sp3Id[4] = "L12";
  static constexpr char rnxId[24] = "CRYOSAT-2";

  /// @ref
  /// https://ids-doris.org/documents/BC/satellites/DORISSatelliteModels.pdf
  static double InitialMass() noexcept {return 724.6e0;} // [kg]

  /// @brief Cryosat-2 initial center of gravity in the satellite fixed frame, 
  /// in cartesian components (x,y,z) [m].
  /// @ref
  /// https://ids-doris.org/documents/BC/satellites/DORISSatelliteModels.pdf
  static Eigen::Matrix<double, 3, 1> initial_center_of_gravity() noexcept {
    constexpr const double com[] = {1.6312e0, 0.0112e0, 0.0137e0};
    return Eigen::Matrix<double, 3, 1>(com);
  }

  /// @brief Get corrections for satellite's mass and center of gravity
  /// coordinates in satellite reference frame.
  ///
  /// The recdords in this file are "corrections" to the mass and [XYZ]
  /// coordinates of the CoG of the satellite in the satellite reference
  /// frame. Hence, the values this function retrievds, should be added to
  /// the initial values of the respective quantities (see InitialMass and
  /// initial_center_of_gravity).
  ///
  /// @param[in] j3mass An IDS/CNES format file with records of
  /// mass & Center of gravity history, see:
  /// https://ids-doris.org/documents/BC/satellites/DORISSatelliteModels.pdf
  /// and ftp://ftp.ids-doris.org/pub/ids/satellites/cs2mass.txt
  /// @param[in] t Current datetime
  /// @param[out] dmdxyz An array holding corrections for datetime t, in
  /// order:
  ///   dmdxyz[0]: Dmass [kg]
  ///   dmdxyz[1,2,3]: Dx, Dy, Dz [m]
  static int corrections(const char *csmass,
                         const dso::datetime<dso::nanoseconds> &t,
                         double *dmdxyz) {
    return get_satellite_corrections(csmass, t, dmdxyz);
  }

  /// @brief Get satellite's mass and center of gravity coordinates in
  /// satellite reference frame.
  static int mass_cog(const dso::datetime<dso::nanoseconds> &t,
                      const char *csmass, double &mass,
                      Eigen::Matrix<double, 3, 1> &cog_xyz) noexcept {

    double cor[4];
    if (corrections(csmass, t, cor))
      return 1;

    mass = InitialMass() + cor[0]; // [kg]

    cog_xyz = initial_center_of_gravity(); // [m]
    cog_xyz(0) += cor[1];
    cog_xyz(1) += cor[2];
    cog_xyz(2) += cor[3];

    return 0;
  }

  /// @brief DORIS phase center for  the 2GHz carrier in satellite reference
  /// frame; cartesian (x,y,z) in [m]
  /// @ref
  /// https://ids-doris.org/documents/BC/satellites/DORISSatelliteModels.pdf
  static int
  pco(Eigen::Matrix<double, 3, 1> &l1_pco,
      Eigen::Matrix<double, 3, 1> &l2_pco) noexcept {
    constexpr const double l1pco[] = {1.848e0, -0.200e0, -0.751e0};
    constexpr const double l2pco[] = {1.832e0, -0.200e0, -0.598e0};
    l1_pco = Eigen::Matrix<double, 3, 1>(l1pco);
    l2_pco = Eigen::Matrix<double, 3, 1>(l2pco);
    return 0;
  }

}; // SatelliteInfo<SATELLITE::Cryosat2>
}// dso

#endif