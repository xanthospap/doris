#include "tides.hpp"
#include "geodesy/geodesy.hpp"
#include <cmath>

namespace {
/// @see IERS 2010 Conventions, Sec. 7.1.4
/// @brief First part of computations, that evaluate m1 and m2 parameters,
///        independent of station coordinates
void pole_tide_part1(double tfyears, double xp, double yp, double &m1,
                     double &m2) noexcept {
  // compute secular pole at t (IERS2010, Sec 7.1.4, Eq. 21) in milliarcseconds
  // transform to arcseconds
  const double xs = (55e0 + 1677e-3 * (tfyears - 2e3)) * 1e-3;
  const double ys = (3205e-1 + 3460e-3 * (tfyears - 2e3)) * 1e-3;

  // assuming xp, yp in arcseconds, apply Eq. 25
  m1 = xp - xs;
  m2 = -(yp - ys);
}

/// @see IERS 2010 Conventions, Sec. 7.1.4
/// @brief Second part of computations, that take m1 and m2 parameters as
///        input, and compute station deiplacements
Eigen::Matrix<double, 3, 1>
pole_tide_part2(double m1, double m2,
                const Eigen::Matrix<double, 3, 1> &r) noexcept {
  // TODO is the following reasoning correct (for the colatitude)?
  // we need the co-latitude, hence we need to convert the ellipsoidal
  // latitude to geocentric and then take the complementary angle
  const double theta =
      dso::DPI / 2e0 -
      dso::geocentric_latitude<dso::ellipsoid::grs80>(r(1), r(2));

  // trigonometric numbers
  const double sl = std::sin(r(0));         // sin(longitude)
  const double cl = std::cos(r(0));         // cos(longitude)
  const double ct = std::cos(theta);        // cos(colatitude)
  const double s2t = std::sin(1e0 * theta); // sin(2*colatitude)
  const double c2t = 2e0 * ct * ct - 1e0;   // cos(2*colatitude)

  // apply Eq. 26 to find the displacements in millimeters
  const double Stheta = -9e0 * c2t * (m1 * cl + m2 * sl);
  const double Slambda = 9e0 * ct * (m1 * sl - m2 * cl);
  const double Sr = -33e0 * s2t * (m1 * cl + m2 * sl);

  return Eigen::Matrix<double, 3, 1>(Slambda, Stheta, Sr);
}
}// unnamed namespace

/// @brief Compute Pole Tides according to IERS 2010, Sec. 7.1.4
/// @param t  Date (IERS does not identify the time-scale)
/// @param xp Polar motion in x [arcseconds]
/// @param yp Polar motion in y [arcseconds]
/// @param sites Geodetic coordinates aka longitude, latitude height (in
///           this order), in units of [rad], [rad], [m]
/// @param Senu The deformation vector as (Slambda, Stheta, Sr) aka in east, 
///           north and up directions, in [mm]. If the site has cartesian 
///           coordinates (x,y,z), then the changes in them due to the pole 
///           tide are: (dx,dy,dz) = R^(T) * (Stheta, Slambda, Sr), with
///             | cosθcosλ  cosθsinλ  -sinθ |
///         R = |  -sinθ      cosλ      0   |
///             | sinθcosλ  sinθsinλ   cosθ |
/// @return Always 0
int dso::pole_tide(const dso::datetime<dso::nanoseconds> &t, double xp,
                   double yp,
                   const std::vector<Eigen::Matrix<double, 3, 1>> &sites,
                   std::vector<Eigen::Matrix<double, 3, 1>> &Senu) noexcept {
  if (Senu.capacity() < sites.size())
    Senu.reserve(sites.size());

        Senu.clear();

  // t to fractional years
  const double tfyears = t.as_fractional_years();

  // first part of computations (compute m1 and m2)
  double m1, m2;
  pole_tide_part1(tfyears, xp, yp, m1, m2);

  // for every station, compute the corrections
  for (const auto &site : sites)
    Senu.emplace_back(pole_tide_part2(m1, m2, site));

  // all done
  return 0;
}

/// @param[in] tfyears Date in years of 365.25 days
/// @param[in] xp Polar motion in x [arcseconds]
/// @param[in] yp Polar motion in y [arcseconds]
/// @param[in] r  Geodetic coordinates aka longitude, latitude height (in
///              this order), in units of [rad], [rad], [m]
/// @return The deformation vector as (Slambda, Stheta, Sr) aka in east, north
///         and up directions, in [mm]. If the site has cartesian coordinates
///         (x,y,z), then the changes in them due to the pole tide are:
///         (dx,dy,dz) = R^(T) * (Stheta, Slambda, Sr), with
///             | cosθcosλ  cosθsinλ  -sinθ |
///         R = |  -sinθ      cosλ      0   |
///             | sinθcosλ  sinθsinλ   cosθ |
/// @see IERS 2010 Conventions, Sec. 7.1.4
Eigen::Matrix<double, 3, 1>
dso::pole_tide(double tfyears, double xp, double yp,
               const Eigen::Matrix<double, 3, 1> &r) noexcept {
  // compute parameters m1 and m2
  double m1, m2;
  pole_tide_part1(tfyears, xp, yp, m1, m2);

  // compute the corrections ...
  return pole_tide_part2(m1, m2, r);
}