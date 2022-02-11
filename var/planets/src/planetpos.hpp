#ifndef __PLANET_LOWP_POSITION_CALCULATOR_HPP__
#define __PLANET_LOWP_POSITION_CALCULATOR_HPP__

#include "datetime/dtcalendar.hpp"

namespace dso{

/// @brief Sun's geocentric position using a low precision analytical series
/// @param[in] tt_jc Terestrial time as Julian cent. since J2000
/// @param[out] rsun Components of the Solar position vector [m] with respect
///             to the mean equator and equinox of J2000 (EME2000, ICRF)
/// @see Montenbruck, Gill, Satellite Orbits, Models Methods Applications;
///      Chapter 3.3.2
int sun_vector(double tt_mjd, double *rsun) noexcept;

#if __cplusplus >= 202002L
template <typename T>
requires(T::is_of_sec_type)
#else
template <class T,
          typename = std::enable_if_t<(T::is_of_sec_type)>>
#endif
int sun_vector(const dso::datetime<T> &t, double *rsun) noexcept {
    return sun_vector(t.jcenturies_sinceJ2000(), rsun);
}

/// @brief Moon's geocentric position using a low precision analytical series
/// @param[in] tt_jc Terestrial time as Julian cent. since J2000
/// @param[out] rmon Components of the Lunar position vector [m] with respect
///             to the mean equator and equinox of J2000 (EME2000, ICRF)
/// @see Montenbruck, Gill, Satellite Orbits, Models Methods Applications;
///      Chapter 3.3.2
/// @note that the accuracy here is several arcminutes for lunar longitude and 
///       latitude and  about 500 km fot the lunar distance.
int moon_vector(double tt_mjd, double *rsun) noexcept;

#if __cplusplus >= 202002L
template <typename T>
requires(T::is_of_sec_type)
#else
template <class T,
          typename = std::enable_if_t<(T::is_of_sec_type)>>
#endif
int moon_vector(const dso::datetime<T> &t, double *rmoon) noexcept {
    return moon_vector(t.jcenturies_sinceJ2000(), rmoon);
}

}// dso

#endif
