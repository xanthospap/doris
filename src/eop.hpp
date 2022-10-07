#ifndef __DSO__IERS_BULLTEIN_PARSERS_HPP__
#define __DSO__IERS_BULLTEIN_PARSERS_HPP__

#include "datetime/dtcalendar.hpp"

namespace dso {

struct EopRecord {
  double mjd, xp, yp, ut1, dx, dy;
};

//int interpolate_eop(double fmjd_utc, const double *mjd, const double *xpa,
//                const double *ypa, const double *ut1a, int sz, double &xp,
//                double &yp, double &dut1) noexcept;

/// @brief A struct to hold EOP information
/// @tparam Capacity The max size of the member arrays. Not to be confused
///         with the actual size of the arrays, which is sz. Elements 
///         with indexes in Capacity-sz are invalid.
struct EopLookUpTable {
  ///< actual size of arrays (<= Capacity)
  int sz{0};
  double *mem_arena{nullptr};
  ///< arrays of EOP values extracted from C04
  double *mjd{nullptr}, // [UTC]
  *xpa {nullptr},  // [mas]
  *ypa {nullptr},  // [mas]
  *ut1a{nullptr},  // [msec]
  *dxa{nullptr}, // [mas]
  *dya{nullptr}; // [mas]

  EopLookUpTable(int sz_=10) noexcept :
    sz{sz_},
    mem_arena(new double[sz_ * 6]),
    mjd(mem_arena),
    xpa(mem_arena + sz_),
    ypa(mem_arena + 2*sz_),
    ut1a(mem_arena + 3*sz_),
    dxa(mem_arena + 4*sz_),
    dya(mem_arena + 5*sz_)
  {}

  void resize(int sz_) noexcept {
    if (sz_ < sz) {
      sz = sz_;
      return;
    }
    delete[] mem_arena;
    mem_arena=new double[sz_ * 6];
    mjd =mem_arena;
    xpa =mem_arena + sz_;
    ypa =mem_arena + 2*sz_;
    ut1a=mem_arena + 3*sz_;
    dxa=mem_arena + 4*sz_;
    dya=mem_arena + 5*sz_;
  }

  ~EopLookUpTable() noexcept {
    if (mem_arena) delete[] mem_arena;
    sz = 0;
  }

/// @brief Interpolate and correct to get the xPole, Ypole and DUT1 values at 
///   a given date
///
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
/// @note Asserts instance sz is at least 4 and units are [mas] and [msec]. 
///       That is INTERP uses a window of 4 data points.
///
/// @param[in] fmjd_utc Interpolation epoch; must be at the samte time-scale 
///            as the mjd input array
/// @param[out] xp Computed xPole at input fmjd_utc [mas]
/// @param[out] yp Computed yPole at input fmjd_utc [mas]
/// @param[out] dut1 Computed UT1-UTC at input fmjd_utc [msec]
/// @return Anything other than 0 is an error
int interpolate(double fmjd_utc, double &xp, double &yp,
                  double &dut1) const noexcept;
int interpolate(double fmjd_utc, double &dx, double &dy) const noexcept;
int interpolate(double fmjd_utc, EopRecord &eopr) const noexcept {
  int error = interpolate(fmjd_utc, eopr.xp, eopr.yp, eopr.ut1);
  error += interpolate(fmjd_utc, eopr.dx, eopr.dy);
  return error;
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
  /// @param[out] eoptable An EopLookUpTable. If needed, it will be resized to 
  ///             the size requested (aka end-start) and will hold the 
  ///             following:
  ///             mjd   Array of resolved MJDs 
  ///             xpa   Array of x pole (EOP) in milliarcsecond [mas]
  ///             ypa   Array of y pole (EOP) in milliarcsecond [mas]
  ///             ut1a  Array of UT1-UTC values in millisecond [ms]
  /// @return Anything other than 0 denotes an error
  int parse(dso::modified_julian_day start, dso::modified_julian_day end,
            EopLookUpTable &eoptable) noexcept;
}; // EopFile

} // dso

#endif
