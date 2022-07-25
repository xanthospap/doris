#ifndef __DSO_DORIS_ATMOSPHERE_MODELS_HPP__
#define __DSO_DORIS_ATMOSPHERE_MODELS_HPP__

#include <cstdint>
#include "atmosphere/nrlmsise00.hpp"

namespace dso {

namespace air_density_models {

namespace exponential {
  /// @brief Calculate Atmospheric Density based on Exponential Model
  /// This simple, static model assumes the density of the atmosphere decays 
  /// exponentially with increasing altitude. It also assumes a spherically 
  /// symmetrical distribution of particles, in which the density, r, varies 
  /// exponentially.
  /// Although a very simple approach, this method yields moderate results for 
  /// general studies.
  /// @param[in] sat_altitude_km Satellite altitude in [km] (note that altitude
  ///            is found by substracting the satellite's radius from the 
  ///            Earth's radius)
  /// @return Atmospheric density at given altitude in [km/m^3]
  /// @see Fundamentals of Astrodynamics, Vallado, Chapter 8.6.2
  double density(double sat_altitude_km) noexcept;
}// namespace exponential

// namespace nrlmsise00 {
/// @brief NRLMSISE-00 Model, 2001
/// The NRLMSIS-00 empirical atmosphere model was developed by Mike
/// Picone, Alan Hedin, and Doug Drob based on the MSISE90 model.
/// The MSISE90 model describes the neutral temperature and densities in
/// Earth's atmosphere from ground to thermospheric heights. Below 72.5 km
/// the model is primarily based on the MAP Handbook (Labitzke et al.,
/// 1985) tabulation of zonal average temperature and pressure by Barnett
/// and Corney, which was also used for the CIRA-86. Below 20 km these
/// data were supplemented with averages from the National Meteorological
/// Center (NMC). In addition, pitot tube, falling sphere, and grenade
/// sounder rocket measurements from 1947 to 1972 were taken into
/// consideration. Above 72.5 km MSISE-90 is essentially a revised MSIS-86
/// model taking into account data derived from space shuttle flights and
/// newer incoherent scatter results. For someone interested only in the
/// thermosphere (above 120 km), the author recommends the MSIS-86
/// model. MSISE is also not the model of preference for specialized
/// tropospheric work. It is rather for studies that reach across several
/// atmospheric boundaries.
/// (quoted from http://nssdc.gsfc.nasa.gov/space/model/atmos/nrlmsise00.html)
/// @see https://ccmc.gsfc.nasa.gov/modelweb/atmos/nrlmsise00.html
/// @param[in]  switches See note below (array of index 24)
/// @param[in]	doy    day of year 
/// @param[in]	sec    fseconds in day (UT) 
/// @param[in]	alt    altitude in kilometers 
/// @param[in]	glat   geodetic latitude 
/// @param[in]	glong  geodetic longitude 
/// @param[in]	lst    local apparent solar time (hours), see note below 
/// @param[in]	f107A  81 day average of F10.7 flux (centered on doy) 
/// @param[in]	f107   daily F10.7 flux for previous day 
/// @param[in]	magnetic_index magnetic index(daily)
/// @param[in]  magnetic_array Array containing the following magnetic values
///             as described below (see also the note on switches; if 
///             switches[9]!=-1, this could be a null pointer)
/// @param[out] outd containing the following values:
///            d[0] - He number density(cm-3)
///            d[1] - O number density(cm-3)
///            d[2] - N2 number density(cm-3)
///            d[3] - O2 number density(cm-3)
///            d[4] - AR number density(cm-3)                       
///            d[5] - total mass density(gm/cm3) [includes d[8] in td7d]
///            d[6] - H number density(cm-3)
///            d[7] - N number density(cm-3)
///            d[8] - Anomalous oxygen number density(cm-3)
/// @param[out] outt containing the following values:
///            t[0] - exospheric temperature
///            t[1] - temperature at alt
/// @note The magnetic_array is an array of size 7, containing the following
///       values:
///  0 : daily AP
///  1 : 3 hr AP index for current time
///  2 : 3 hr AP index for 3 hrs before current time
///  3 : 3 hr AP index for 6 hrs before current time
///  4 : 3 hr AP index for 9 hrs before current time
///  5 : Average of eight 3 hr AP indicies from 12 to 33 hrs 
///          prior to current time
///  6 : Average of eight 3 hr AP indicies from 36 to 57 hrs 
///          prior to current time
/// 
/// @note O, H, and N are set to zero below 72.5 km
/// t[0], Exospheric temperature, is set to global average for
/// altitudes below 120 km. The 120 km gradient is left at global
/// average value for altitudes below 72 km.
/// 
/// d[5], TOTAL MASS DENSITY, is NOT the same for subroutines GTD7 
/// and GTD7D
/// 
/// SUBROUTINE GTD7 -- d[5] is the sum of the mass densities of the
/// species labeled by indices 0-4 and 6-7 in output variable d.
/// This includes He, O, N2, O2, Ar, H, and N but does NOT include
/// anomalous oxygen (species index 8).
/// 
/// SUBROUTINE GTD7D -- d[5] is the "effective total mass density
/// for drag" and is the sum of the mass densities of all species
/// in this model, INCLUDING anomalous oxygen.
///
/// @note UT, Local Time, and Longitude are used independently in the
/// model and are not of equal importance for every situation.  
/// For the most physically realistic calculation these three
/// variables should be consistent (lst=sec/3600 + g_long/15).
/// The Equation of Time departures from the above formula
/// for apparent local time can be included if available but
/// are of minor importance.
/// 
/// f107 and f107A values used to generate the model correspond
/// to the 10.7 cm radio flux at the actual distance of the Earth
/// from the Sun rather than the radio flux at 1 AU. The following
/// site provides both classes of values:
/// ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
/// 
/// f107, f107A, and ap effects are neither large nor well
/// established below 80 km and these parameters should be set to
/// 150., 150., and 4. respectively.
///
/// @note Switches: to turn on and off particular variations use these switches.
/// 0 is off, 1 is on, and 2 is main effects off but cross terms on.
/// 
/// Standard values are 0 for switch 0 and 1 for switches 1 to 23. The 
/// array "switches" needs to be set accordingly by the calling program. 
/// The arrays sw and swc are set internally.
/// 
/// switches[i]:
///  i - explanation
/// -----------------
///  0 - output in meters and kilograms instead of centimeters and grams
///  1 - F10.7 effect on mean
///  2 - time independent
///  3 - symmetrical annual
///  4 - symmetrical semiannual
///  5 - asymmetrical annual
///  6 - asymmetrical semiannual
///  7 - diurnal
///  8 - semidiurnal
///  9 - daily ap [when this is set to -1 (!) the pointer amgnetic_array must
///      hold the respective values; else it is not used
/// 10 - all UT/long effects
/// 11 - longitudinal
/// 12 - UT and mixed UT/long
/// 13 - mixed AP/UT/LONG
/// 14 - terdiurnal
/// 15 - departures from diffusive equilibrium
/// 16 - all TINF var
/// 17 - all TLB var
/// 18 - all TN1 var
/// 19 - all S var
/// 20 - all TN2 var
/// 21 - all NLB var
/// 22 - all TN3 var
/// 23 - turbo scale height var
/*void gtd7(
    const int *switches, int doy, double fsec, double glat,
    double glon, double lst, double f107, double f107A, double alt,
    double magnetic_index, const double *magnetic_array, 
    double *outd, double *outt) noexcept;*/
//} // namespace nrlmsise00

namespace msis86 {
  /// @brief MSIS-86/cira 1986 neutral thermosphere model
  /// @ref Hedin, A. E., MSIS-86 thermospheric model, Journal of Geophysical 
  ///      Research, vol. 92, pp. 4649–4662, 1987. doi:10.1029/JA092iA05p04649.
  /// @param[in] day day number of the year.
  /// @param[in] sec ut [sec]
  /// @param[in] alt  altitude [km] (greater than 85 km)
  /// @param[in] glat geodetic latitude [deg]
  /// @param[in] glong geodetic longitude [deg]
  /// @param[in] stl local apparent solar time [hrs]
  /// @param[in] f107a 3 month average of f10.7 flux
  /// @param[in] f107 daily f10.7 flux for previous day
  /// @param[in] ap magnetic index (daily) or when sw(9)=-1. an array 
  ///               containing:
  ///          (0) daily ap
  ///          (1) 3 hr ap index for current time
  ///          (2) 3 hr ap index for 3 hrs before current time
  ///          (3) 3 hr ap index for 6 hrs before current time
  ///          (4) 3 hr ap index for 9 hrs before current time
  ///          (5) average of eight 3 hr ap indicies from 12 to 33 hrs prior
  ///              to current time
  ///          (6) average of eight 3 hr ap indicies from 36 to 59 hrs prior
  ///              to current time
  /// @param[out] d
  ///   d(0) - He number density [cm^3]
  ///   d(1) - O number density [cm^3]
  ///   d(2) - N2 number density [cm^3]
  ///   d(3) - O2 number density [cm^3]
  ///   d(4) - Ar number density [cm^3]
  ///   d(5) - total mass density (gm/cm3)
  ///   d(6) - H number density [cm^3]
  ///   d(7) - N number density[cm^3]
  /// @param[out] t
  ///   t(0) - exospheric temperature
  ///   t(1) - temperature at alt
  ///          (2) to turn on and off particular variations call tselec(sw)
  /// @note switches is a 25 element array containing 0. for off, 1. for on, or 
  /// 2. for main effects off but cross terms on for the following variations
  ///     1 - f10.7 effect on mean  2 - time independent
  ///     3 - symmetrical annual    4 - symmetrical semiannual
  ///     5 - asymmetrical annual   6 - asymmetrical semiannual
  ///     7 - diurnal               8 - semidiurnal
  ///     9 - daily ap             10 - all ut/long effects
  ///    11 - longitudinal         12 - ut and mixed ut/long
  ///    13 - mixed ap/ut/long     14 - terdiurnal
  ///    15 - departures from diffusive equilibrium
  ///    16 - all tinf var         17 - all tlb var
  ///    18 - all t0 var           19 - all s var
  ///    20 - all z0 var           21 - all nlb var
  ///    22 - all tr12 var         23 - turbo scale height var
  ///
  ///    - Fortran implementation can be found here: 
  ///      https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/msis86/
  ///    - Matlab implementation 
  ///      https://www.researchgate.net/publication/337085334_MSIS-86_Atmosphere_Model_MATLAB_code
  int msis86(int doy, double sec, double alt, double glat, double glong,
           double stl, double f107a, double f107, const int *switches,
           const double *ap, double *d, double *t) noexcept;
} // msis86
} // air_density_models
} // dso

#endif
