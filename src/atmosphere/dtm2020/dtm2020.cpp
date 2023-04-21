#include "atmosphere/dtm2020/dtm2020.hpp"
#include "geodesy/geoconst.hpp"
#include <cstdio>
#include <stdexcept>

double dso::details::DIN__day_of_year(const dso::TwoPartDate &t) noexcept {
  const auto ydoy =
      (dso::modified_julian_day{(long)t.normalized().big()}).to_ydoy();
  double doy = (double)(ydoy.__doy.as_underlying_type());
  return doy + t.normalized().small();
}

/* TODO is this correct ? */
double dso::details::DIN__local_hours_radians(const dso::TwoPartDate &t,
                                              double lon) noexcept {
  constexpr const double _1RadHours = 24e0 / dso::D2PI;
  constexpr const double _SecInHour = 60 * 60e0;
  /* difference from UTC at Greenwich in fractional hours */
  const double dhours =
      lon * _1RadHours + std::copysign(t.normalized().small() * _SecInHour, lon);
  /* add difference in Grewnwich time if lon > 0 subtract the delta hours  */
  return (t.normalized().small() * _SecInHour) +
         (std::copysign(dhours, lon) * (-1e0));
}

dso::Dtm2020::Dtm2020(const char *fn) : in{}, out{} {
  if (map_model_coeffs(fn)) {
    fprintf(stderr,
            "[ERROR] Failed to parse DTM2020 data file %s. Aborting "
            "construction (traceback: %s)\n",
            fn, __func__);
    throw std::runtime_error("FAILED DTM2020 CTOR\n");
  }
}
