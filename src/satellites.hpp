#ifndef __DSO_DORIS_SATELLITES_HPP__
#define __DSO_DORIS_SATELLITES_HPP__

#include "datetime/dtcalendar.hpp"
#include "eigen3/Eigen/Eigen"

namespace dso {

/* @brief Get corrections for satellite mass and center of gravity
 * coordinates in satellite reference frame.
 *
 * The recdords in this file are "corrections" to the mass and [XYZ]
 * coordinates of the CoG of the satellite in the satellite reference
 * frame. Hence, the values this function retrievds, should be added to
 * the initial values of the respective quantities (see InitialMass and
 * initial_center_of_gravity).
 *
 * @param[in] j3mass An IDS/CNES format file with records of
 * mass & Center of gravity history for Jason-3 see:
 * https://ids-doris.org/documents/BC/satellites/DORISSatelliteModels.pdf
 * (e.g. ftp://ftp.ids-doris.org/pub/ids/satellites/ja3mass.txt for Jason-3)
 * @param[in] t Current datetime
 * @param[out] dmdxyz An array holding corrections for datetime t, in
 * order:
 *   dmdxyz[0]: Dmass [kg]
 *   dmdxyz[1,2,3]: Dx, Dy, Dz [m]
 */
int get_satellite_corrections(const char *ids_mass,
                              const dso::datetime<dso::nanoseconds> &t,
                              double *dmdxyz) noexcept;

struct MacroModelComponent {
  /* Surface [m^2] */
  double m_surf;
  double area() const noexcept {return m_surf;}
  /* Normal in satellite fixed ref. frame [x,y,z] */
  double m_normal[3];
  Eigen::Matrix<double, 3, 1> normal() const noexcept {
    return Eigen::Matrix<double, 3, 1>(m_normal);
  }
  void set_normal(const Eigen::Matrix<double, 3, 1> &n) noexcept {
    m_normal[0] = n(0);
    m_normal[1] = n(1);
    m_normal[2] = n(2);
  }
  /* spec, diff, abs */
  double m_optical_properties[3];
  double optical_specular() const noexcept {return m_optical_properties[0];}
  double optical_diffusion() const noexcept {return m_optical_properties[1];}
  double optical_absorption() const noexcept {return m_optical_properties[2];}
  /* spec, diff, abs */
  double m_infrared_properties[3];
};

enum class SATELLITE : char { Jason3, Cryosat2 }; // SATELLITE
template <SATELLITE Sat> struct SatelliteInfo {};
template <SATELLITE Sat> struct MacroModel {};

} // namespace dso

#endif
