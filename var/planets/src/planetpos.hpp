#ifndef __PLANET_LOWP_POSITION_CALCULATOR_HPP__
#define __PLANET_LOWP_POSITION_CALCULATOR_HPP__

namespace dso{
/// @brief Sun's geocentric position using a low precision analytical series
/// @param[in] tt_mjd Terestrial time as Modified Julian Date
/// @param[out] rsun Components of the Solar position vector [m] with respect
///             to the mean equator and equinox of J2000 (EME2000, ICRF)
/// @see Montenbruck, Gill, Satellite Orbits, Models Methods Applications;
///      Chapter 3.3.2
int sun_vector(double tt_mjd, double *rsun) noexcept;

/// @brief Moon's geocentric position using a low precision analytical series
/// @param[in] tt_mjd Terestrial time as Modified Julian Date
/// @param[out] rmon Components of the Lunar position vector [m] with respect
///             to the mean equator and equinox of J2000 (EME2000, ICRF)
/// @see Montenbruck, Gill, Satellite Orbits, Models Methods Applications;
///      Chapter 3.3.2
int moon_vector(double tt_mjd, double *rsun) noexcept;
}// dso

#endif
