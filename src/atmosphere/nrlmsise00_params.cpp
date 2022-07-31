#include "nrlmsise00.hpp"

dso::nrlmsise00::detail::InParamsCore::InParamsCore(dso::modified_julian_day mjd,
                                      double secday) noexcept {
  const dso::datetime<dso::milliseconds> t(
      mjd, dso::milliseconds(
               static_cast<dso::milliseconds::underlying_type>(secday * 1e3)));
  const auto ydoy = t.as_ydoy();
  year = ydoy.__year.as_underlying_type();
  doy = ydoy.__doy.as_underlying_type();
  sec = t.sec().to_fractional_seconds();
}
