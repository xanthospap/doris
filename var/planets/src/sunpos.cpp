#include "planetpos.hpp"
#include "geodesy/units.hpp"
#include "celgeo/iau.hpp"

using std::sin;
using std::cos;

int dso::sun_vector(double t, double *rsun) noexcept {
  /*const double t =
      (tt_mjd - 51544e5) / 36525e0; */

    constexpr double f1 = dso::deg2rad<double>(357.5256e0);
    constexpr double f2 = dso::deg2rad<double>(35999.049e0);
    const double M = f1 + f2 * t;

    constexpr double omega_sum = dso::deg2rad<double>(282.94e0);
    const double sinM = sin(M);
    const double sin2M = sin(2e0*M);
    const double cosM = cos(M);
    const double cos2M = cos(2e0*M);
    
    // Sun's ecliptic logitude (rad)
    const double lambda = omega_sum + M + 6892e0 * sinM + 72e0 * sin2M;
    // distance (m)
    const double r = 149.619e0 - 2.499e0 * cosM - 0.021e0 * cos2M;

    constexpr double ecliptic_obliquity = dso::deg2rad<double>(23.43929111);
    constexpr double sinEps = sin(ecliptic_obliquity);
    constexpr double cosEps = cos(ecliptic_obliquity);
    const double sinl = sin(lambda);
    const double cosl = cos(lambda);
    
    rsun[0] = r * cosl;
    rsun[1] = r * sinl * cosEps;
    rsun[2] = r * sinl * sinEps;

    return 0;
}

