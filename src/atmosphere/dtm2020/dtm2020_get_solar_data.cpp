#include "atmosphere/dtm2020/dtm2020.hpp"
#include <cassert>

using dso::utils::celestrak::details::MissingSwData;
constexpr const double RECOMPUTE_AFTER_MIN = 3;

int dso::Dtm2020::get_solar_data(const dso::TwoPartDate &tutcin,
                                 bool use_adjusted_values) noexcept {
  const dso::TwoPartDate tutc = tutcin.normalized();

  /* do not redo this shit if we have alredy computed it up to Xmin berfore */
  static dso::TwoPartDate LastUsedUTC = dso::TwoPartDate(36934, 0e0);
  if (double sdt = std::abs(
          tutc.diff<dso::DateTimeDifferenceType::FractionalSeconds>(LastUsedUTC));
      sdt < RECOMPUTE_AFTER_MIN * 3600e0) {
    return 0;
  }
  
  /* t - 24hours */
  const dso::TwoPartDate t24(tutc._big - 1, tutc._small);

  /* instantaneous flux at (t - 24hr); interpolate between two days */
  double flux24 = MissingSwData;
  {
    const dso::TwoPartDate t1mjd(t24._big, 0);
    const dso::TwoPartDate t2mjd(t24._big + 1, 0);
    /* indexes in array */
    auto it = flux_data.at(t1mjd);
    if ((it == flux_data.data_end()) || (it + 1 == flux_data.data_end())) {
      fprintf(stderr,
              "[ERROR] Failed to interpolate F107 data! Date is out of range "
              "(traceback: %s)\n",
              __func__);
      return 1;
    }
    /* linear interpolation */
    const dso::TwoPartDate t1((double)it->mjd_.as_underlying_type(), 0);
    const double y1 = it->f107Obs * (!use_adjusted_values) +
                      it->f107Adj * (use_adjusted_values);
    const dso::TwoPartDate t2((double)(it + 1)->mjd_.as_underlying_type(), 0);
    const double y2 = (it + 1)->f107Obs * (!use_adjusted_values) +
                      (it + 1)->f107Adj * (use_adjusted_values);
    assert(t1<t24 && t2>=t24);
    const double dt1 = t24.diff<dso::DateTimeDifferenceType::FractionalDays>(t1);
    const double dt2 = t2.diff<dso::DateTimeDifferenceType::FractionalDays>(t1);
    flux24 = y1 + dt1 * (y2 - y1) / dt2;
    if (y1 == MissingSwData || y2 == MissingSwData) {
      fprintf(stderr, "[ERROR] Missing Solar flux data (traceback: %s)\n",
              __func__);
      return 1;
    }
  }

  /* mean flux of last 81 days at t */
  const double flux81 =
      flux_data.at(tutc)->f107ObsL81 * (!use_adjusted_values) +
      flux_data.at(tutc)->f107AdjL81 * (use_adjusted_values);
  if (flux81 == MissingSwData) {
    fprintf(stderr, "[ERROR] Missing Solar flux data (traceback: %s)\n",
            __func__);
    return 1;
  }

  /* akp(1)= kp delayed by 3 hours */
  double akp1;
  {
    const dso::TwoPartDate t3h(
        dso::TwoPartDate(tutc._big, tutc._small - (3e0 / 24)).normalized());
    auto it = flux_data.at(t3h);
    if (it == flux_data.data_end()) {
      fprintf(stderr,
              "[ERROR] Failed to interpolate F107 data! Date is out of range "
              "(traceback: %s)\n",
              __func__);
      return 1;
    }
    akp1 = it->KpIndexes[static_cast<int>((t3h._small * 24) / 3)];
  }

  /* akp(3)=mean of last 24 hours */
  double akp24=MissingSwData;
  {
    auto it = flux_data.at(t24);
    if (it == flux_data.data_end()) {
      fprintf(stderr,
              "[ERROR] Failed to interpolate F107 data! Date is out of range "
              "(traceback: %s)\n",
              __func__);
      return 1;
    }
    int error = 0;
    int start_index = static_cast<int>((t24._small * 24) / 3);
    for (int indexes_read = 0; indexes_read < 8; indexes_read++) {
      error += (it->KpIndexes[start_index] == MissingSwData);
      akp24 += it->KpIndexes[start_index++];
      if (start_index == 8) {
        start_index = 0;
        ++it;
      }
    }
    if (error) {
      fprintf(stderr,
              "[ERROR] Failed to interpolate F107 data! Date is out of range "
              "(traceback: %s)\n",
              __func__);
      return 1;
    }
    akp24 /= 8e0;
  }

  /* set last computed datetime */
  LastUsedUTC = tutc;

  /* assign to arrays */
  mf[0] = flux24;
  mf[1] = 0e0;
  mfbar[0] = flux81;
  mfbar[1] = 0e0;
  makp[0] = akp1;
  makp[1] = 0e0;
  makp[2] = akp24;
  makp[3] = 0e0;

  /* all done */
  return 0;
}
