#ifndef __DSO_DORIS_ATMOSPHERE_MODELS_HPP__
#define __DSO_DORIS_ATMOSPHERE_MODELS_HPP__

#include <cstdint>

namespace dso {

namespace air_density_models {

namespace exponential {

}// namespace exponential

namespace nrlmsise00 {
/// @brief NRLMSISE-00 Model, 2001
/// @see https://ccmc.gsfc.nasa.gov/modelweb/atmos/nrlmsise00.html
} // namespace nrlmsise00

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
