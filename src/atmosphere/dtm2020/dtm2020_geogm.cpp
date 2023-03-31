#include "atmosphere/dtm2020/dtm2020.hpp"
#include "geodesy/units.hpp"
#include <cmath>
#include <cstdio>

/* Original FORTRAN routine for this function found at
 * https://github.com/swami-h2020-eu/mcm/blob/main/src/dtm2020/dtm2020_F30_Hp-subr.f90
 *
 * WARNING
 * The original FORTRAN version handles angles (input/output) using dec. 
 * degrees. This version uses radians!
 */

int dso::details::geogm(double xlat, double xlong, double &glat,
                        double &glon) noexcept {

  if (xlong == dso::rad2deg(291e0))
    xlong = dso::deg2rad(291.1e0);
  constexpr const double platr = dso::deg2rad(78.5);
  constexpr const double plongr = dso::deg2rad(291e0);

  const double spl = std::sin(platr);
  const double cpl = std::cos(platr);
  const double slm =
      spl * std::sin(xlat) + cpl * std::cos(xlat) * std::cos(plongr - xlong);
  const double clm = std::sqrt(1. - slm * slm);
  const double phim1 = std::cos(xlat) * std::sin(xlong - plongr) / clm;
  const double phim2 = (spl * slm - std::sin(xlat)) / (cpl * clm);
  glat = std::asin(slm);

  /* to degrees for checks */
  const double phim1_deg = dso::rad2deg(phim1);
  const double phim2_deg = dso::rad2deg(phim2);

  /* please make this brachless, its a mess! TODO */
  if ((phim1_deg >= 0e0) && (phim2_deg >= 0e0)) {
    glon = std::asin(phim1);
    return 0;
  }
  if ((phim1_deg >= 0e0) && (phim2_deg < 0e0)) {
    glon = dso::DPI - std::abs(std::asin(phim1));
    return 0;
  }
  if ((phim1_deg < 0e0) && (phim2_deg < 0e0)) {
    glon = dso::DPI + std::abs(std::asin(phim1));
    return 0;
  }
  if ((phim1_deg < 0e0) && (phim2_deg >= 0e0)) {
    glon = dso::D2PI - std::abs(std::asin(phim1));
    return 0;
  }

    fprintf(stderr,
            "[ERROR] Failed converting geographic to geomagnetic coordinates! "
            "(traceback: %s)\n",
            __func__);
    return 1;
}
