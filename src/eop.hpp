#ifndef __DSO__IERS_BULLTEIN_PARSERS_HPP__
#define __DSO__IERS_BULLTEIN_PARSERS_HPP__

#include "datetime/dtcalendar.hpp"

namespace dso {

/// @brief Given EOP data (of size sz), interpolate and correct to get the
///        xPole, Ypole and DUT1 values at a given date
/// The function will perform the following:
///   1. first use INTERP to interpolate x,y,ut1 values for the given date.
///   INTERP is recommended to interpolate the IERS polar motion and Universal 
///   Time products and account for the semidiurnal/diurnal variations in the 
///   Earths orientation. This procedure makes use of a Lagrangian 
///   interpolation scheme and applies the integral Ray model (71 tidal waves) 
///   and Brzezinski-Mathews-Bretagnon-Capitaine-Bizouard model (10 lunisolar 
///   waves) of the semidiurnal/diurnal variations in the Earth's orientation 
///   as recommended in the IERS 2000 Conventions (McCarthy, 2002). see
///   https://hpiers.obspm.fr/iers/models/interp.readme
///   2. Use ORTHO_EOP to ccount for variations in polar motion (Dx,Dy) 
///   due to ocean-tides
///   3. Use PMSDNUT2 to account for libration effects
/// See IERS Conventions 2010, Petit et al., Chapter 5.5
///
/// @note INTERP (aka iers2010::interp_pole` uses a window of 4 points to
///       interpolate. That means we should have data for at least one day 
///       before fmjd_utc and tow days after it.
///
/// @param[in] fmjd_utc Interpolation epoch; must be at the samte time-scale 
///            as the mjd input array
/// @param[in] mjd MJD dates of input data (size sz)
/// @param[in] xpa xpole data at input mjd (size sz) [mas]
/// @param[in] ypa xpole data at input mjd (size sz) [mas]
/// @param[in] dut1 UT1-UTC data at input mjd (size sz) [mas]
/// @param[in] sz Size of input arrays
/// @param[out] xp Computed xPole at input fmjd_utc [mas]
/// @param[out] yp Computed yPole at input fmjd_utc [mas]
/// @param[out] dut1 Computed UT1-UTC at input fmjd_utc [msec]
/// @return Anything other than 0 is an error
int interpolate_eop(double fmjd_utc, const double *mjd, const double *xpa,
                const double *ypa, const double *ut1a, int sz, double &xp,
                double &yp, double &dut1) noexcept;

/// @brief A struct to hold EOP information
/// @tparam Capacity The max size of the member arrays. Not to be confused
///         with the actual size of the arrays, which is sz. Elements 
///         with indexes in Capacity-sz are invalid.
template<int Capacity>
struct EopLookUpTable {
  ///< actual size of arrays (<= Capacity)
  int sz;
  ///< arrays of EOP values extracted from C04
  double mjd[Capacity], // [UTC]
  xpa[Capacity],  // [mas]
  ypa[Capacity],  // [mas]
  ut1a[Capacity]; // [msec]

  /// @brief Interpolates and applies ocean tide and libration corrections to
  ///        EOP data for given date.
  /// @see dso::interpolate_eop
  int interpolate(double fmjd_utc, double &xp, double &yp,
                  double &dut1) const noexcept {
    return interpolate_eop(fmjd_utc, mjd, xpa, ypa, ut1a, sz, xp, yp, dut1);
  }
};// EopLookUpTable

/// @brief EopFile is a (dead simple) file EOP information, just like
///        IERS Bulletin B/C04 files. It's actually a translation of such 
///        files, but holding only EOP data for dates of choice.
class EopFile {
  ///< filename
  char filename[256];

public:
  /// @brief Constructor from filename
  EopFile(const char *fn);
  /// @brief Copy not allowed
  EopFile(const EopFile &other) = delete;
  /// @brief Move allowed
  EopFile(EopFile &&other) noexcept;
  /// @brief Assignment not allowed
  EopFile &operator=(const EopFile &other) = delete;
  /// @brief Move assignment allowed
  EopFile &operator=(EopFile &&other) noexcept;
  /// @brief Destructor
  ~EopFile() noexcept {};

  /// @brief Extract data EOP from an EopFile for given dates
  ///        The function will extract EOP for the time interval: [start,end) 
  ///        off of this instance. 
  ///        The data will be stored in the passed in arrays, which must be
  ///        of size >= end - start
  ///
  /// @warning Asserts that data in the file are in chronological order.
  ///
  /// @param[in]  start MJD for start date (included)
  /// @param[in]  end   MJD for end date (not included)
  /// @param[out] mjd   Array of resolved MJDs 
  /// @param[out] xpa   Array of x pole (EOP) in milliarcsecond [mas]
  /// @param[out] ypa   Array of y pole (EOP) in milliarcsecond [mas]
  /// @param[out] ut1a  Array of UT1-UTC values in millisecond [ms]
  /// @return Anything other than 0 denotes an error
  int parse(dso::modified_julian_day start, dso::modified_julian_day end,
            double *mjd, double *xpa, double *ypa, double *ut1a) noexcept;

  /// @brief Extract data EOP from an EopFile for given dates
  ///        The function will extract EOP for the time interval: [start,end) 
  ///        off of this instance. 
  ///        The data will be stored in the passed in arrays, which must be
  ///        of size >= end - start
  template <int Capacity>
  int make_lookup_table(dso::modified_julian_day start,
                        dso::modified_julian_day end,
                        EopLookUpTable<Capacity> &eop_lut) noexcept {
    int sz = end.as_underlying_type() - start.as_underlying_type();
    if (sz < Capacity)
      return 1;
    const int error = parse(start, end, eop_lut.mjd, eop_lut.xpa, eop_lut.ypa, eop_lut.ut1a);
    eop_lut.sz = (error==0) ? sz : 0;
    return error;
  }
}; // EopFile

} // dso

#endif
