#ifndef __DSO__IERS_BULLTEIN_PARSERS_HPP__
#define __DSO__IERS_BULLTEIN_PARSERS_HPP__

#include "datetime/dtcalendar.hpp"
#include <cstring>

namespace dso {

/// @brief A simple record of EOP values
struct EopRecord {
  double mjd, xp, yp, dut, lod, dx, dy, omega;
};

/// @brief A class to hold EOP information ordered by Mjd
///        The elements of the EOP/MJD arrays are always stored in
///        chronological order. Note that notmally, (i.e. if the instance has
///        been constructed by a call to dso::parse_iers_C04) the mjd array is
///        in TT (whereas the IERS-distributed C04 files contain date
///        information in UTC).
class EopLookUpTable {
private:
  ///< Actual size of arrays aka number of epochs
  int sz{0};
  ///< Capacity of memmory, not the same as the actual number of epochs stored!
  int capacity{0};
  ///< Memmory pool
  double *mem_arena{nullptr};
  ///< arrays of EOP values extracted from C04
  double *mjda{nullptr}, ///< [TT], start index = 0
      *xpa{nullptr},     ///< [", arcsec], start index = 1
      *ypa{nullptr},     ///< [", arcsec], start index = 2
      *ut1a{nullptr},    ///< [sec], start index = 3
      *dxa{nullptr},     ///< [", arcsec], start index = 4
      *dya{nullptr},     ///< [", arcsec], start index = 5
      *loda{nullptr},    ///< [sec/day], start index = 6
      *omegaa{nullptr};  ///< [?], start index = 7

public:
  int size() const noexcept { return sz; }

  double *mjd(int i = 0) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return mjda + i;
  }
  double *xp(int i = 0) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return xpa + i;
  }
  double *yp(int i = 0) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return ypa + i;
  }
  double *dut(int i = 0) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return ut1a + i;
  }
  double *dx(int i = 0) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return dxa + i;
  }
  double *dy(int i = 0) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return dya + i;
  }
  double *lod(int i = 0) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return loda + i;
  }
  double *omega(int i = 0) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return omegaa + i;
  }
  const double *mjd(int i = 0) const noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return mjda + i;
  }
  const double *xp(int i = 0) const noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return xpa + i;
  }
  const double *yp(int i = 0) const noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return ypa + i;
  }
  const double *dut(int i = 0) const noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return ut1a + i;
  }
  const double *dx(int i = 0) const noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return dxa + i;
  }
  const double *dy(int i = 0) const noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return dya + i;
  }
  const double *lod(int i = 0) const noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return loda + i;
  }
  const double *omega(int i = 0) const noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < sz);
#endif
    return omegaa + i;
  }

  /// @brief Default constructor (uses a capacity of 10 elements)
  EopLookUpTable(int _capacity = 10) noexcept
      : sz{_capacity}, capacity{_capacity},
        mem_arena(new double[_capacity * 8]), mjda(mem_arena),
        xpa(mem_arena + _capacity), ypa(mem_arena + 2 * _capacity),
        ut1a(mem_arena + 3 * _capacity), dxa(mem_arena + 4 * _capacity),
        dya(mem_arena + 5 * _capacity), loda(mem_arena + 6 * _capacity),
        omegaa(mem_arena + 7 * _capacity) {}

  /// @brief Resize; alocation/dealocation depends on capacity
  void resize(int sz_) noexcept {
    // no need to dealocate/alocate; just change the effective size
    if (sz_ < capacity) {
      sz = sz_;
      return;
    } else {
      // requested size (sz_) > current capacity. need to allocate memory!
      delete[] mem_arena;
      sz = sz_;
      capacity = sz_;
      mem_arena = new double[capacity * 8];
      mjda = mem_arena;
      xpa = mem_arena + capacity;
      ypa = mem_arena + 2 * capacity;
      ut1a = mem_arena + 3 * capacity;
      dxa = mem_arena + 4 * capacity;
      dya = mem_arena + 5 * capacity;
      loda = mem_arena + 6 * capacity;
      omegaa = mem_arena + 7 * capacity;
    }
  }

  /// @brief destructor
  ~EopLookUpTable() noexcept {
    if (mem_arena)
      delete[] mem_arena;
    sz = capacity = 0;
  }

  /// @brief Move constructor
  EopLookUpTable(EopLookUpTable &&eopt) noexcept
      : sz{eopt.sz}, capacity{eopt.sz}, mem_arena{eopt.mem_arena},
        mjda{eopt.mjda}, xpa{eopt.xpa}, ypa{eopt.ypa}, ut1a{eopt.ut1a},
        dxa{eopt.dxa}, dya{eopt.dya}, loda{eopt.loda}, omegaa{eopt.omegaa} {
    eopt.sz = eopt.capacity = 0;
    eopt.mem_arena = nullptr;
  }

  /// @brief Assignment operator
  EopLookUpTable &operator=(const EopLookUpTable &eopt) noexcept {
    if (this != &eopt) {
      if (sz != eopt.sz) {
        this->resize(eopt.sz);
      }
      mjda = (double *)std::memcpy(mem_arena, eopt.mem_arena, sz);
      xpa = (double *)std::memcpy(mem_arena + 1 * capacity,
                                  eopt.mem_arena + 1 * eopt.capacity, sz);
      ypa = (double *)std::memcpy(mem_arena + 2 * capacity,
                                  eopt.mem_arena + 2 * eopt.capacity, sz);
      ut1a = (double *)std::memcpy(mem_arena + 3 * capacity,
                                   eopt.mem_arena + 3 * eopt.capacity, sz);
      dxa = (double *)std::memcpy(mem_arena + 4 * capacity,
                                  eopt.mem_arena + 4 * eopt.capacity, sz);
      dya = (double *)std::memcpy(mem_arena + 5 * capacity,
                                  eopt.mem_arena + 5 * eopt.capacity, sz);
      loda = (double *)std::memcpy(mem_arena + 6 * capacity,
                                   eopt.mem_arena + 6 * eopt.capacity, sz);
      omegaa = (double *)std::memcpy(mem_arena + 7 * capacity,
                                     eopt.mem_arena + 7 * eopt.capacity, sz);
    }
    return *this;
  }

  /// @brief Move assignment operator
  EopLookUpTable &operator=(EopLookUpTable &&eopt) noexcept {
    if (this != &eopt) {
      // just move elements ...
      sz = eopt.sz;
      capacity = eopt.sz;
      mem_arena = eopt.mem_arena;
      mjda = eopt.mjda;
      xpa = eopt.xpa;
      ypa = eopt.ypa;
      ut1a = eopt.ut1a;
      dxa = eopt.dxa;
      dya = eopt.dya;
      loda = eopt.loda;
      omegaa = eopt.omegaa;
      // leave move-from instance blank
      eopt.sz = 0;
      eopt.capacity = 0;
      eopt.mem_arena = nullptr;
    }
    return *this;
  }

  /// @brief Copy constructor. The created copy will have a capacity and size
  ///        (both) equal to eopt's size.
  EopLookUpTable(const EopLookUpTable &eopt) noexcept
      : sz{eopt.sz}, capacity{eopt.sz}, mem_arena(new double[capacity * 8]),
        mjda((double *)std::memcpy(mem_arena, eopt.mem_arena, sz)),
        xpa((double *)std::memcpy(mem_arena + 1 * capacity,
                                  eopt.mem_arena + 1 * eopt.capacity, sz)),
        ypa((double *)std::memcpy(mem_arena + 2 * capacity,
                                  eopt.mem_arena + 2 * eopt.capacity, sz)),
        ut1a((double *)std::memcpy(mem_arena + 3 * capacity,
                                   eopt.mem_arena + 3 * eopt.capacity, sz)),
        dxa((double *)std::memcpy(mem_arena + 4 * capacity,
                                  eopt.mem_arena + 4 * eopt.capacity, sz)),
        dya((double *)std::memcpy(mem_arena + 5 * capacity,
                                  eopt.mem_arena + 5 * eopt.capacity, sz)),
        loda((double *)std::memcpy(mem_arena + 6 * capacity,
                                   eopt.mem_arena + 6 * eopt.capacity, sz)),
        omegaa((double *)std::memcpy(mem_arena + 7 * capacity,
                                     eopt.mem_arena + 7 * eopt.capacity, sz)) {}

  /// @brief Compute the effect of zonal Earth tides on the rotation of the
  ///        Earth (using iers2010::rg_zont2 on UT1-UTC, LOD and
  ///        Omega) and subtract it from the respective elements in the
  ///        instance (hence this value will change the ut1a, loda and omegaa
  ///        arrays. This will result on the so-called "regularized" EOPs.
  void __regularize() noexcept;

  /// @brief Interpolate and correct to get EOP/ERP parameters at given date
  ///
  /// See IERS Conventions 2010, Petit et al., Chapter 5.5
  ///
  /// @param[in] fmjd_tt Interpolation epoch given as mjd [TT] (assuming that
  ///            the calling instance's mjd array is in TT)
  /// @param[out] eopr EOP/ERP parameters at requested epoch, as computed by
  ///            Lagrangian interpolation and the removal of ocen-tide effects
  ///            (iers2010::interp::pmut1_oceans) and libration
  ///            (iers2010::interp::pm_gravi) effects.
  /// @return Anything other than 0 is an error
  /// TODO Must handle leap seconds, see Bradley et al
  int interpolate(double fmjd_tt, EopRecord &eopr,
                  int order = 3) const noexcept;
  int interpolate2(double fmjd_tt, EopRecord &eopr,
                  int order = 3) const noexcept;

  /// @brief Interpolate EOP/ERPs at time ffmjd_tt using Lgrangian
  /// interpolation
  ///        of given order
  /// @param fmjd_tt Point to interpolate at [MJD] TT
  /// @param eopr Interpolation results stored in an EopRecord instance;
  /// note
  ///             that omega (Earth rotation rate) values are NOT
  ///             used/interpolated
  /// @param order Order of the Lagrangian interpolation (that is the
  /// window, aka
  ///             used points is order+1). This parameter should be an odd
  ///             integer
  /// @return Anything other than 0 denotes an error
  int interpolate_lagrange(double fmjd_tt, EopRecord &eopr,
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
/// @note  EOP files contain UTC time-stamps. In this function, the time-tags 
///        (mjd array) will be converted to TT before appended to the mjd 
///        array, hence the mjd array in the eoptable instance will be n 
///        MJD/TT.
/// @return Anything other than 0 denotes an error
int parse_iers_C04(const char *c04fn, dso::modified_julian_day start,
                   dso::modified_julian_day end,
                   EopLookUpTable &eoptable) noexcept;
} // namespace dso

#endif
