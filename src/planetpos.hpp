#ifndef __PLANET_LOWP_POSITION_CALCULATOR_HPP__
#define __PLANET_LOWP_POSITION_CALCULATOR_HPP__

#include "base_error.hpp"
#include "cppspice/SpiceUsr.h"
#include "datetime/dtcalendar.hpp"
#include "eigen3/Eigen/Eigen"

/* @todo I have not tested run speeds, but i should make the approximate
 * functions as efficient as possible cause then, what is the point? Users
 * could extract sun/moon coordinates from planetary ephemeris (cspice)
 */
namespace dso {

/* Enum class to denote kinda planets, for easy interaction with CSPICE */
enum class Planet : char {
  EARTH,
  MOON,
  SUN,
  MERCURY,
  VENUS,
  MARS,
  JUPITER,
  SATURN
};

namespace cspice {

/* @brief Transform a Planet enum to its NAIF Id */
dso::iStatus planet_to_naif_id(dso::Planet p, int &id) noexcept;

/* @brief Load a cspice kernel (of any type)
 * @param[in] kernel The name of the cspice kernel, see
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/kernel.html#Kernel%20Types
 * @note This is a wrapper function around furnsh_c; see
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/furnsh_c.html
 */
inline int load_kernel(const char *kernel) noexcept {
  furnsh_c(kernel);
  return 0;
}

/* @brief Load a given ascii/LSK cspice kernel if not already loaded.
 * @see
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/kernel.html#Kernel%20Types
 */
int load_if_unloaded_lsk(const char *lsk_kernel) noexcept;

/* @brief Load a given binary/SPK cspice kernel if not already loaded.
 * @see
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/kernel.html#Kernel%20Types
 */
int load_if_unloaded_spk(const char *spk_kernel) noexcept;

/* @brief Load a given binary/PCK cspice kernel if not already loaded.
 * @see
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/kernel.html#Kernel%20Types
 */
int load_if_unloaded_pck(const char *pck_kernel) noexcept;

/* @brief Compute Ephemeris Time from TT passed as Julian Date via cspice
 * @note The function expects that:
 *       * An lsk (aka leap-second) kernel is already loaded, so that we
 *       can perform datetime computations
 */
inline double mjdtt2et(const dso::TwoPartDate &tt) noexcept {
  /* get ephemeris (ET) time from TT; assumes an LSK kernel is loaded ... */
  return unitim_c(tt.jd(), "JDTDT", "ET");
}

/* @brief Position of a target body relative to an observing body.
 * The reference frame for the returned vector is J2000 and the (returned)
 * values have units [km]. Note that no aberration corrections are applied.
 * @note The function expects that:
 *       * An spk kernel is already loaded so that we can compute the
 *         target/oberver positions
 * For the target/oberver ID's, see
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html
 * Note that according to CSPICE documentation, J200 frame is actually ICRF,
 * see
 * https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/04_concepts.pdf
 *
 * @param[in] et The epoch as ephemeris time (see mjdtt2et)
 * @param[in] target_id An int representing the target body, using the
 *            interbal scpice IDs (see planet_to_naif_id)
 * @param[in] observer_id An int representing the observing body, using the
 *            interbal scpice IDs (planet_to_naif_id)
 * @param[out] pos A vector of size >= 3, where the X/Y/Z components of the
 *            vector are stored; units [km]
 */
inline int j2planet_pos_from(double et, int target_id, int observer_id,
                             double *pos) noexcept {
  double dummy;
  /* get the position of the planet */
  spkezp_c(target_id, et, "J2000", "NONE", observer_id, pos, &dummy);
  return 0;
}

} /* namespace cspice */

/* @brief get planet position w.r.t Earth in GCRF/J2000
 * @pamra[in] p A Planet enum, denoting the target planet
 * @param[in] mjd_tt The date of request as an MJD in TT
 * @param[out] pos The poistion of the (target) planet at the date of request
 *            in [m] at the GCRF/J2000 frame
 */
iStatus planet_pos(Planet p, const TwoPartDate &mjd_tt,
                   Eigen::Matrix<double, 3, 1> &pos) noexcept;

/* @brief Get the Sun's and Moon's standard gravitational parameter (μ) off
 *        from a CSPICE/NAIF PCK Kernel
 * @param[in] pck_kernel A PCK kernel filename. If the kernel is already
 *         loaded, we will get the values without reloading it. If the
 *         filename is NULL, then we will try retrieving the constants
 *         assuming a PCK kernel is already loaded.
 * @param[out] GMSun Sun's gravitational constant according to the kernel in
 *         [km^3 / sec^2]
 * @param[out] GMMoon Moon's gravitational constant according to the kernel in
 *         [km^3 / sec^2]
 * @param[in] use_si Transform Sun's and Moon's gravitational constant to SI
 *         units, i.e. to m^3/sec^2 (instead of km^3/sec^2)
 * @return Anything other than zero denotes an error.
 */
int get_sun_moon_GM(const char *pck_kernel, double &GMSun,
                    double &GMMoon, int use_si=false) noexcept;

} /* namespace dso */

#endif
