#ifndef __DSO__IERS_BULLTEIN_PARSERS_HPP__
#define __DSO__IERS_BULLTEIN_PARSERS_HPP__

#include "datetime/dtcalendar.hpp"
#include <datetime/dtfund.hpp>
#include <fstream>
#include <limits>

namespace dso {

/// Bulletin-B/C04 contain up to one month of final data (aka before the
/// day of publication) and up to one month of prelimenery (prediction) data
/// after the publication date. We can use these values as max sizes when 
/// parsing data from these files.

namespace bulletin_details {
constexpr const double BULLETIN_MISSING_VALUE = 999.999e0;
constexpr const int FILE_IS_AHEAD = std::numeric_limits<int>::min();
constexpr const int FILE_IS_PRIOR = std::numeric_limits<int>::min() + 1;
} // namespace bulletin_details

template<int N>
struct EopLookUpTable {
  ///< actual size of arrays (<= N)
  int sz;
  ///< arrays of EOP values extracted from C04
  double mjd[N], xpa[N], ypa[N], ut1a[N];
};// EopLookUpTable

struct [[deprecated]] IersBulletinB_Section1Block {
  long mjd;
  double x, y, dut1, dX, dY; //[mas], [mas], [ms], [mas], [mas]
  double xerr, yerr, dut1err, dXerr, dYerr;
  double dUt1rUt = bulletin_details::BULLETIN_MISSING_VALUE;
  char type; // 'F' for final, 'P' for preliminery
};

class EopFile {
  char filename[256];
  EopFile(const char *fn);
  EopFile(const EopFile &other) = delete;
  EopFile(EopFile &&other) noexcept;
  EopFile &operator=(const EopFile &other) = delete;
  EopFile &operator=(EopFile &&other) noexcept;
  ~EopFile() noexcept;

  //
  int parse(dso::modified_julian_day start, dso::modified_julian_day end,
            double *mjd, double *xpa, double *ypa, double *ut1a, int &sz) noexcept;
}; // EopFile

class [[deprecated]] IersBulletinB {
private:
  char filename[256];
  std::ifstream stream;
  bool has_dUt1UtcR;

public:
  IersBulletinB(const char *fn);
  IersBulletinB(const IersBulletinB &other) = delete;
  IersBulletinB(IersBulletinB &&other) noexcept;
  IersBulletinB &operator=(const IersBulletinB &other) = delete;
  IersBulletinB &operator=(IersBulletinB &&other) noexcept;
  ~IersBulletinB() noexcept;

  /// @brief Return true if the file has UT1R-UTC values
  bool has_dUt1UtcR_info() const noexcept { return this->has_dUt1UtcR; }

  /// @brief Parse an IERS Bulletin B (aka IersBulletinB isntance) Section 1
  /// records
  /// @param[out] block An array of IersBulletinB_Section1Block where the
  ///             parsed Section 1 record lines are stored at (must be large 
  ///             enough to hold as many blocks as needed). The function will 
  ///             return the number of blocks filled within the array.
  /// @param[in]  include_preliminary Signals if preliminary values are to
  ///             be collected or not.
  /// @return An integer; if < 0 an error has ocured. If >=0 the number of
  ///         lines parsed; hence, the number of blocks filled in the input
  ///         block array.
  int parse_section1(IersBulletinB_Section1Block *block,
                     bool include_preliminary = true) noexcept;

  /// @brief Search through an IERS Bulletin B (aka IersBulletinB instance)
  /// Section 1 block, to find records for the given date. Bulletin B files
  /// only record integral days, hence only integer parts of MJDs are compared
  /// (input datetime t's fractional part is of no significance).
  /// @param[in]  mjd The MJD which we are searching for; only the integral
  ///             part is neede, hence an int is provided.
  /// @param[out] block At output, holds the Section 1 record for the given date
  /// @param[in]  include_preliminary Signals if preliminary values are to be
  ///             collected or not.
  /// @return Anything other than 0, denotes an error. If the returned value is
  /// < 0, then the given date was not found. If it is  > 0, a parsing error
  /// occured.
  /// @warning Be careful on how you use the extracted values. According to
  ///          IERS Conventions, e.g. x/y pole coordinates need to be corrected
  ///          for a number of effects before used in Terrestrial/Celestial
  ///          transformations.
  int get_section1_at(int mjd, IersBulletinB_Section1Block &block,
                      bool include_preliminary = true) noexcept;

#if __cplusplus >= 202002L
  template <gconcepts::is_sec_dt S>
#else
  template <typename S, typename = std::enable_if_t<S::is_of_sec_type>>
#endif
  inline int get_section1_at(const datetime<S> &t,
                             IersBulletinB_Section1Block &block,
                             bool include_preliminary = true) noexcept {
    int mjd = t.mjs().as_underlying_type();
    return get_section1_at(mjd, block, include_preliminary);
  }

}; // IersBulletinB

[[deprecated]]
int download_iers_bulletinb_for(long mjd, char *downloaded_fn = nullptr,
                                const char *dir = nullptr) noexcept;

} // namespace dso

#endif
