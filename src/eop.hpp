#ifndef __DSO__IERS_BULLTEIN_PARSERS_HPP__
#define __DSO__IERS_BULLTEIN_PARSERS_HPP__

#include "datetime/dtcalendar.hpp"
#include <cstring>
#include <vector>

namespace dso {

/// @brief A simple record of EOP values
struct EopRecord {
  dso::TwoPartDate mjd;
  double xp, yp, dut, lod, dx, dy, omega;
};

/// @brief A class to hold EOP information ordered by Mjd
///        The elements of the EOP/MJD arrays are always stored in
///        chronological order. Note that the mjd array can contain dates in
///        either TT or UTC depending on how the IERS-distributed C04 file has
///        been parsed (see parse_iers_C04)
struct EopLookUpTable {
  std::vector<dso::TwoPartDate> t;
  std::vector<double> xp, yp, dut1, dX, dY, lod;

  EopLookUpTable() noexcept { reserve(8); }

  void push_back(const EopRecord &r) noexcept {
    t.push_back(r.mjd);
    xp.push_back(r.xp);
    yp.push_back(r.yp);
    dut1.push_back(r.dut);
    dX.push_back(r.dx);
    dY.push_back(r.dy);
    lod.push_back(r.lod);
  }

  int size() const noexcept { return t.size(); }

  void clear() noexcept {
    t.clear();
    xp.clear();
    yp.clear();
    dut1.clear();
    dX.clear();
    dY.clear();
    lod.clear();
  }

  void reserve(int sz) noexcept {
    t.reserve(sz);
    xp.reserve(sz);
    yp.reserve(sz);
    dut1.reserve(sz);
    dX.reserve(sz);
    dY.reserve(sz);
    lod.reserve(sz);
  }

  /// @brief Interpolate and correct to get EOP/ERP parameters at given date
  ///
  /// See IERS Conventions 2010, Petit et al., Chapter 5.5
  ///
  /// @param[in] mjd Interpolation epoch given as mjd [UTC] (assuming that
  ///            the calling instance's mjd array is in UTC)
  /// @param[out] eopr EOP/ERP parameters at requested epoch, as computed by
  ///            Lagrangian interpolation and the removal of ocen-tide effects
  ///            (iers2010::interp::pmut1_oceans) and libration
  ///            (iers2010::interp::pm_gravi) effects.
  /// @return Anything other than 0 is an error
  /// TODO Must handle leap seconds, see Bradley et al
  int interpolate(const dso::TwoPartDate &mjd, EopRecord &eopr,
                  int order = 3) const noexcept;

  /// @brief Interpolate EOP/ERPs at time ffmjd_tt using Lgrangian
  /// interpolation
  ///        of given order
  /// @param mjd Point to interpolate at [MJD] UTC
  /// @param eopr Interpolation results stored in an EopRecord instance;
  /// note
  ///             that omega (Earth rotation rate) values are NOT
  ///             used/interpolated
  /// @param order Order of the Lagrangian interpolation (that is the
  /// window, aka
  ///             used points is order+1). This parameter should be an odd
  ///             integer
  /// @return Anything other than 0 denotes an error
  int interpolate_lagrange(const dso::TwoPartDate &mjd, EopRecord &eopr,
                           int order = 3) const noexcept;
}; // EopLookUpTable

/// @brief Extract data EOP from an EopFile for given dates
///        The function will extract EOP for the time interval: [start,end)
///        off of this instance.
///        The data will be stored in the passed in arrays, which must be
///        of size >= end - start
///
/// @warning Asserts that data in the file are in chronological order.
///
/// @param[in]  start MJD for start date (included) [UTC]
/// @param[in]  end   MJD for end date (not included) [UTC]
/// @param[out] eoptable An EopLookUpTable. If needed, it will be resized to
///             the size requested (aka end-start) and will hold the
///             following:
///             mjd   Array of resolved MJDs [TT]
///             xpa   Array of x pole (EOP) in arcseconds [arcsec]
///             ypa   Array of y pole (EOP) in arcseconds [arcsec]
///             ut1a  Array of UT1-UTC values in seconds [sec]
///             loda  Array of LOD values in seconds/day [sec]
/// @param[utc2tt] EOP files contain UTC time-stamps. In set to true, the
///             time-tags (mjd array) will be converted to TT before appended
///             to the mjd array, hence the mjd array in the eoptable instance
///             will be in MJD/TT.
/// @return Anything other than 0 denotes an error
int parse_iers_C04(const char *c04fn, dso::modified_julian_day start,
                   dso::modified_julian_day end, EopLookUpTable &eoptable,
                   bool utc2tt = false) noexcept;
} // namespace dso

#endif
