#ifndef __PLANET_LOWP_POSITION_CALCULATOR_HPP__
#define __PLANET_LOWP_POSITION_CALCULATOR_HPP__

#include "datetime/dtcalendar.hpp"
#include "cppspice/SpiceUsr.h"

/// @todo I have not tested run speeds, but i should make the approximate
/// functions as efficient as possible cause then, what is the point? Users
/// could extract sun/moon coordinates from planetary ephemeris (cspice)

namespace dso {

namespace cspice {
  /// @brief Load a cspice kernel (of any type)
  /// @param[in] kernel The name of the cspice kernel, see
  /// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/kernel.html#Kernel%20Types
  /// @note This is a wrapper function around furnsh_c; see
  /// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/furnsh_c.html
  inline
  int load_kernel(const char *kernel) noexcept {
    furnsh_c(kernel);
    return 0;
  }

  /// @brief Load a given ascii/LSK cspice kernel if not already loaded.
  /// @see https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/kernel.html#Kernel%20Types
  int load_if_unloaded_lsk(const char *lsk_kernel) noexcept;
  
  /// @brief Load a given binary/SPK cspice kernel if not already loaded.
  /// @see https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/kernel.html#Kernel%20Types
  int load_if_unloaded_spk(const char *spk_kernel) noexcept;
  
  /// @brief Load a given binary/PCK cspice kernel if not already loaded.
  /// @see https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/kernel.html#Kernel%20Types
  int load_if_unloaded_pck(const char *pck_kernel) noexcept;

  /// @brief Compute Ephemeris Time from TT passed as Julian Date via cspice
  /// @note The function expects that:
  ///       * An lsk (aka leap-second) kernel is already loaded, so that we
  ///         can perform datetime computations
  inline
  double jd2et(double jd_tt) noexcept {
    // get ephemeris (ET) time from TT; assumes an LSK kernel is loaded ...
    return unitim_c(jd_tt, "JDTDT", "ET");
  }

  /// @brief Position of a target body relative to an observing body.
  /// The reference frame for the returned vector is J2000 and the (returned)
  /// values have units [km]. Note that no aberration corrections are applied.
  /// @note The function expects that:
  ///       * An spk kernel is already loaded so that we can compute the
  ///         target/oberver positions
  /// For the target/oberver ID's, see 
  /// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
  ///
  /// @param[in] et The epoch as ephemeris time (see jd2et)
  /// @param[in] target_id An int representing the target body, using the 
  ///            interbal scpice IDs
  /// @param[in] observer_id An int representing the observing body, using the 
  ///            interbal scpice IDs
  /// @param[out] pos A vector of size >= 3, where the X/Y/Z components of the
  ///            vector are stored; units [km]
  inline
  int j2planet_pos_from(double et, int target_id, int observer_id,
                        double *pos) noexcept {
    double dummy;
    // get the position of the planet
    spkezp_c(target_id, et, "J2000", "NONE", observer_id, pos, &dummy);
    return 0;
  }
}// spice

/// @brief Get the Sun's and Moon's standard gravitational parameter (μ) off
///        from a CSPICE/NAIF PCK Kernel
/// @param[in] pck_kernel A PCK kernel filename. If the kernel is already
///         loaded, we will get the values without reloading it. If the
///         filename is NULL, then we will try retrieving the constants
///         assuming a PCK kernel is already loaded.
/// @param[out] GMSun Sun's gravitational constant according to the kernel in
///         [km^3 / sec^2]
/// @param[out] GMMoon Moon's gravitational constant according to the kernel in
///         [km^3 / sec^2]
/// @return Anything other than zero denotes an error.
int get_sun_moon_GM(const char *pck_kernel, double &GMSun,
                 double &GMMoon) noexcept;

/// @brief Sun's geocentric position using a low precision analytical series
/// @param[in] tt_jdc Terestrial time as Modified Julian cent. since J2000
/// @param[out] rsun Components of the Solar position vector [km] with respect
///             to the mean equator and equinox of J2000 (EME2000, ICRF)
/// @see Montenbruck, Gill, Satellite Orbits, Models Methods Applications;
///      Chapter 3.3.2
int sun_vector_montenbruck(double tt_jdc, double *rsun) noexcept;


/// @brief Sun's geocentric position using a low precision analytical series
/// @param[in] tt_jdc Terestrial time as Modified Julian cent. since J2000
/// @param[out] rsun Components of the Solar position vector [km] with respect
///             to the mean equator and equinox of J2000 (EME2000, ICRF)
/// @see Vallado, Fundamentals of Astrodynamics and Applications, Fourth 
///             Edition, Chapter 5.1.1
/// @note This is the same algorithm as also described in Orbital Mechanics for
/// Engineering Students, H. D. Curtis, (Fourth Edition), Chapter 10.10
int sun_vector_vallado(double tt_jdc, double *rsun) noexcept;
#ifdef DEBUG
int sun_vector_approx21(double tt_jdc, double *rsun) noexcept;
#endif


/// @brief Sun's geocentric position using a low precision analytical series
/// @param[in] tt_jdc Terestrial time as Modified Julian cent. since J2000
/// @param[out] rsun Components of the Solar position vector [km] with respect
///             to the mean equator and equinox of J2000 (EME2000, ICRF)
/// @see Curtis, Orbital Mechanics for Engineering Students, Fourth Edition,
///             Chapter 12.9
int sun_vector_curtis(double n, double *rsun) noexcept;

/// @brief Sun's geocentric position using planetary ephemeris via cspice
/// The reference frame for the returned vector is J2000 and the (returned)
/// values have units [km]. Note that no aberration corrections are applied.
/// @note The function expects that:
///       * An lsk (aka leap-second) kernel is already loaded, so that we
///         can perform datetime computations, and
///       * An spk kernel is already loaded so that we can compute the
///         target/oberver positions
/// @param[in] t An datetime in TT
/// @param[out] pos A vector of size >= 3, where the X/Y/Z components of the
///            vector are stored; units [km]
#if __cplusplus >= 202002L
template <typename T>
requires(T::is_of_sec_type)
#else
template <class T,
          typename = std::enable_if_t<(T::is_of_sec_type)>>
#endif
int sun_vector_cspice(const dso::datetime<T> &t, double *rsun) noexcept {
  const double jd = t.as_jd(); // date as jd (TT)
  return cspice::j2planet_pos_from(cspice::jd2et(jd), 10, 399, rsun);
}

#if __cplusplus >= 202002L
template <typename T>
requires(T::is_of_sec_type)
#else
template <class T,
          typename = std::enable_if_t<(T::is_of_sec_type)>>
#endif
int sun_vector_montenbruck(const dso::datetime<T> &t, double *rsun) noexcept {
    return sun_vector_montenbruck(t.jcenturies_sinceJ2000(), rsun);
}

#if __cplusplus >= 202002L
template <typename T>
requires(T::is_of_sec_type)
#else
template <class T,
          typename = std::enable_if_t<(T::is_of_sec_type)>>
#endif
int sun_vector_vallado(const dso::datetime<T> &t, double *rsun) noexcept {
    return sun_vector_vallado(t.jcenturies_sinceJ2000(), rsun);
}

#if __cplusplus >= 202002L
template <typename T>
requires(T::is_of_sec_type)
#else
template <class T,
          typename = std::enable_if_t<(T::is_of_sec_type)>>
#endif
int sun_vector_curtis(const dso::datetime<T> &t, double *rsun) noexcept {
    return sun_vector_curtis(t.as_jd() - dso::j2000_jd, rsun);
}

/// @brief Moon's geocentric position using a low precision analytical series
/// @param[in] tt_jc Terestrial time as Julian cent. since J2000
/// @param[out] rmon Components of the Lunar position vector [km] with respect
///             to the mean equator and equinox of J2000 (EME2000, ICRF)
/// @see Montenbruck, Gill, Satellite Orbits, Models Methods Applications;
///      Chapter 3.3.2
/// @note that the accuracy here is several arcminutes for lunar longitude and 
///       latitude and  about 500 km fot the lunar distance.
int moon_vector_approx(double tt_mjd, double *rsun) noexcept;

/// @brief Moon's equatoria geocentric position using an approximate formula
/// @param[in] tdb_jc TDB time as Julian cent. since J2000
/// @param[out] rmon Components of the Lunar position vector [km]
/// @see Vallado, Fundamentals of Astrodynamics and Applications, Fourth 
///             Edition, Chapter 5.2
int moon_vector_vallado(double tt_mjd, double *rsun) noexcept;

#if __cplusplus >= 202002L
template <typename T>
requires(T::is_of_sec_type)
#else
template <class T,
          typename = std::enable_if_t<(T::is_of_sec_type)>>
#endif
int moon_vector_approx(const dso::datetime<T> &t, double *rmoon) noexcept {
    return moon_vector_approx(t.jcenturies_sinceJ2000(), rmoon);
}
#if __cplusplus >= 202002L
template <typename T>
requires(T::is_of_sec_type)
#else
template <class T,
          typename = std::enable_if_t<(T::is_of_sec_type)>>
#endif
int moon_vector_vallado(const dso::datetime<T> &t, double *rmoon) noexcept {
    return moon_vector_vallado(t.jcenturies_sinceJ2000(), rmoon);
}

/// @brief Moon's geocentric position using planetary ephemeris via cspice
/// The reference frame for the returned vector is J2000 and the (returned)
/// values have units [km]. Note that no aberration corrections are applied.
/// @note The function expects that:
///       * An lsk (aka leap-second) kernel is already loaded, so that we
///         can perform datetime computations, and
///       * An spk kernel is already loaded so that we can compute the
///         target/oberver positions
/// @param[in] t An datetime in TT
/// @param[out] pos A vector of size >= 3, where the X/Y/Z components of the
///            vector are stored; units [km]
#if __cplusplus >= 202002L
template <typename T>
requires(T::is_of_sec_type)
#else
template <class T,
          typename = std::enable_if_t<(T::is_of_sec_type)>>
#endif
int moon_vector_cspice(const dso::datetime<T> &t, double *rmoon) noexcept {
  const double jd = t.as_mjd() + dso::mjd0_jd; // date as jd (TT)
  return cspice::j2planet_pos_from(cspice::jd2et(jd), 301, 399, rmoon);
}

}// dso

#endif
