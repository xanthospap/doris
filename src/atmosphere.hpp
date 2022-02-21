#ifndef __DSO_DORIS_ATMOSPHERE_MODELS_HPP__
#define __DSO_DORIS_ATMOSPHERE_MODELS_HPP__

#include <cstdint>

namespace dso {

namespace air_density_models {
namespace nrlmsise00 {
/// @brief NRLMSISE-00 Model, 2001
/// @see https://ccmc.gsfc.nasa.gov/modelweb/atmos/nrlmsise00.html

/// Switches: to turn on and off particular variations use these switches.
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
///  9 - daily ap [when this is set to -1 (!) the pointer
///                ap_a in struct nrlmsise_input must
///                point to a struct ap_array]
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
struct nrlmsise_flags {
  int8_t switches[24] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
}; // nrlmsise_flags

struct flag_arena {
  double SW[24*2];
};

/// @brief Array containing the following magnetic values:
///  0 : daily AP
///  1 : 3 hr AP index for current time
///  2 : 3 hr AP index for 3 hrs before current time
///  3 : 3 hr AP index for 6 hrs before current time
///  4 : 3 hr AP index for 9 hrs before current time
///  5 : Average of eight 3 hr AP indicies from 12 to 33 hrs
///          prior to current time
///  6 : Average of eight 3 hr AP indicies from 36 to 57 hrs
///          prior to current time
//struct ap_array {
//  double a[7];
//};

/// OUTPUT VARIABLES:
///    d[0] - HE NUMBER DENSITY(CM-3)
///    d[1] - O NUMBER DENSITY(CM-3)
///    d[2] - N2 NUMBER DENSITY(CM-3)
///    d[3] - O2 NUMBER DENSITY(CM-3)
///    d[4] - AR NUMBER DENSITY(CM-3)
///    d[5] - TOTAL MASS DENSITY(GM/CM3) [includes d[8] in td7d]
///    d[6] - H NUMBER DENSITY(CM-3)
///    d[7] - N NUMBER DENSITY(CM-3)
///    d[8] - Anomalous oxygen NUMBER DENSITY(CM-3)
///    t[0] - EXOSPHERIC TEMPERATURE
///    t[1] - TEMPERATURE AT ALT
///
///    O, H, and N are set to zero below 72.5 km
///
///    t[0], Exospheric temperature, is set to global average for
///    altitudes below 120 km. The 120 km gradient is left at global
///    average value for altitudes below 72 km.
///
///    d[5], TOTAL MASS DENSITY, is NOT the same for subroutines GTD7
///    and GTD7D
///
///      SUBROUTINE GTD7 -- d[5] is the sum of the mass densities of the
///      species labeled by indices 0-4 and 6-7 in output variable d.
///      This includes He, O, N2, O2, Ar, H, and N but does NOT include
///      anomalous oxygen (species index 8).
///
///      SUBROUTINE GTD7D -- d[5] is the "effective total mass density
///      for drag" and is the sum of the mass densities of all species
///      in this model, INCLUDING anomalous oxygen.
struct nrlmsise_output {
  double d[9]; ///< densities
  double t[2]; ///< temperatures
};             // nrlmsise_output

struct nrlmsise_input {
  int year;      ///< year, currently ignored
  int doy;       ///< day of year
  double sec;    ///< seconds in day (UT)
  double alt;    ///< altitude in kilometers
  double g_lat;  ///< geodetic latitude
  double g_lon;  ///< geodetic longitude
  double lst;    ///< local apparent solar time (hours), see note below
  double f107A;  ///< 81 day average of F10.7 flux (centered on doy)
  double f107;   ///< daily F10.7 flux for previous day
  double ap;     ///< magnetic index(daily)
  double ap_array[7] = {0e0}; ///< see above
};               // nrlmsise_input

void gts7(const nrlmsise_input &input, const nrlmsise_flags &flags,
          const flag_arena &f_arena, nrlmsise_output &output) noexcept;
//void ghp7(const nrlmsise_input &input, nrlmsise_flags &flags,
//          nrlmsise_output &output, double press) noexcept;
void gtd7d(const nrlmsise_input &input, const nrlmsise_flags &flags,
           flag_arena &f_arena, nrlmsise_output &output) noexcept;
void gtd7(const nrlmsise_input &input, const nrlmsise_flags &flags,
          /*flag_arena &f_arena,*/ nrlmsise_output &output) noexcept;

} // namespace nrlmsise00
} // namespace air_density_models
} // namespace dso

#endif
