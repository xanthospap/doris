#ifndef __DSO_DORIS_SATELLITES_HPP__
#define __DSO_DORIS_SATELLITES_HPP__

#include "datetime/dtcalendar.hpp"

namespace dso {

/// @brief Get corrections for satellite mass and center of gravity
/// coordinates in satellite reference frame.
///
/// The recdords in this file are "corrections" to the mass and [XYZ]
/// coordinates of the CoG of the satellite in the satellite reference
/// frame. Hence, the values this function retrievds, should be added to
/// the initial values of the respective quantities (see InitialMass and
/// initial_center_of_gravity).
///
/// @param[in] j3mass An IDS/CNES format file with records of
/// mass & Center of gravity history for Jason-3 see:
/// https://ids-doris.org/documents/BC/satellites/DORISSatelliteModels.pdf
/// (e.g. ftp://ftp.ids-doris.org/pub/ids/satellites/ja3mass.txt for Jason-3)
/// @param[in] t Current datetime
/// @param[out] dmdxyz An array holding corrections for datetime t, in
/// order:
///   dmdxyz[0]: Dmass [kg]
///   dmdxyz[1,2,3]: Dx, Dy, Dz [m]
int get_satellite_corrections(const char *ids_mass,
                           const dso::datetime<dso::nanoseconds> &t,
                           double *dmdxyz) noexcept;

enum class SATELLITE : char { Jason3, Cryosat2 }; // SATELLITE
template <SATELLITE Sat> struct SatelliteInfo {};

} // dso

#endif