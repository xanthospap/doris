#ifndef __SOLAR_ACTIVITY_DATA_FEED_HPP__
#define __SOLAR_ACTIVITY_DATA_FEED_HPP__

#include "var/celestrak.hpp"
#include <vector>
#include <algorithm>

namespace dso {

class SolarActivityData {
  typedef dso::utils::celestrak::details::CelestTrakSWFlux CelestTrakSWFlux;
  std::vector<CelestTrakSWFlux> vdata;

public:
  /* @brief Extract data from the input CelestTrack SCV file for the interval:
   * [t.MJD-days_before, t.MJD+days_after)
   * @param[in] fnSwCsv CelestTrack SCV data file, see 
   *                     https://celestrak.org/SpaceData/
   * @param[in] t       Central date, only its MJD is considered (UTC)
   * @param[in] days_before Specify starting date for data collection. The 
   *                    starting date is actually a day, specified by its MJD
   *                    as starting_day = t.MJD-days_before. Records for this
   *                    MJD are included in the vdata vector.
   * @param[in] days_after Specify (exclusive) ending date for data collection. 
   *                    The ending date is actually a day, specified by its MJD
   *                    as ending_day = t.MJD+days_after. Records for this MJD 
   *                    are NOT included in the vdata vector.
   * @return Anything other than 0 denotes an error
   */
  int feed(const char *fnSwCsv, const dso::TwoPartDate &tutc, int days_before,
           int days_after) noexcept;

  auto at(const TwoPartDate &tutc) const noexcept {
    const TwoPartDate t(tutc.normalized());
    const dso::modified_julian_day cmjd(t.big());
    return std::find_if(
        vdata.begin(), vdata.end(),
        [=](const CelestTrakSWFlux &flux) { return flux.mjd_ == cmjd; });
  }

  auto data_end() const noexcept {
    return vdata.end();
  }
  auto data_begin() const noexcept {
    return vdata.begin();
  }

  const std::vector<CelestTrakSWFlux> &data() const noexcept {return vdata;}
}; // SolarActivityData
} // namespace dso

#endif
