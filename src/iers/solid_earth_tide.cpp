#include "geodesy/geodesy.hpp"
#include "iers2010/iau.hpp"
#include "iers2010/iersc.hpp"
#include "harmonic_coeffs.hpp"
#include "egravity.hpp"
#include "tides.hpp"

namespace {
/// @param[in] dC   Normalized corrections to C coefficients, in the order:
///            dC = C20,C21,C22,C30,C31,C32,C33,C40,C41,C42,0,0
///                   0   1   2   3   4   5   6   7   8   9
/// @param[in] dS   Normalized corrections to S coefficients, in the order:
///            dS = 0,S21,S22,0,S31,S32,S33,0,S41,S42,0,0
///                 0,  1   2 3   4   5   6 7   8   9
int array2harmonics(const std::array<double, 12> &dC,
                    const std::array<double, 12> &dS,
                    dso::HarmonicCoeffs &cs) noexcept 
{
  cs.clear();
  cs.C(2,0) = dC[0];
  cs.C(3,0) = dC[3];
  cs.C(4,0) = dC[7];
  cs.C(2,1) = dC[1];
  cs.C(3,1) = dC[4];
  cs.C(4,1) = dC[8];
  cs.C(2,2) = dC[2];
  cs.C(3,2) = dC[5];
  cs.C(4,2) = dC[9];
  cs.C(3,3) = dC[6];

  cs.S(2,1) = dS[1];
  cs.S(2,2) = dS[2];
  cs.S(3,1) = dS[4];
  cs.S(3,2) = dS[5];
  cs.S(3,3) = dS[6];
  cs.S(4,1) = dS[8];
  cs.S(4,2) = dS[9];

return 0;
}

void compute_legendre(
    double latitude,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &P) noexcept {
  const double sf = std::sin(latitude);
  const double sf2 = sf*sf;
  const double cf = std::cos(latitude);
  const double cf2 = cf*cf;
  P(0,0) = 0e0;

  P(1,0) = sf;
  P(1,1) = -cf;

  P(2,0) = 0.5e0 * (3e0 * sf2 - 1e0);
  P(2,1) = -3e0 * sf *cf;
  P(2,2) = 3e0 * cf;

  P(3,0) = 0.5e0*sf*(5e0*sf2-3e0);
  P(3,1) = (3e0/2e0)*(1e0-5e0*sf2)*cf;
  P(3,2) = 15e0*sf*cf2;
  P(3,3) = -15e0*cf*cf2;

  P(4,0) = (1e0/8e0)*(35e0*sf2*sf2-30e0*sf2+3e0);
  P(4,1) = (-5e0/2e0)*(7e0*sf2*sf-3e0*sf)*cf;
  P(4,2) = 7.5e0 * (7e0*sf*sf-1e0)*cf2;
  P(4,3) = -105e0*sf*cf2*cf;
  P(4,4) = 105e0*cf2*cf2;

  return;
}
}// unnamed namespace

/// @brief Compute corrections to normalized C and S gravitational
///        coefficients, using the model(s) described in IERS2010 standards.
///
///        The corrections are returned in the dC and dS (ouput) arrays.
///
/// @param[in] t_tt  datetime in TT
/// @param[in] ut1_mjd Corresponding datetime in UT1 time scale, as MJD
/// @param[in] rmoon ECEF vector of Moon, as (X,Y,Z) in [m]
/// @param[in] rsun E CEF vector of Sun, as (X,Y,Z) in [m]
/// @param[out] dC   Normalized corrections to C coefficients, in the order:
///             dC = C20,C21,C22,C30,C31,C32,C33,C40,C41,C42,0,0
/// @param[out] dS   Normalized corrections to S coefficients, in the order:
///             dS = 0,S21,S22,S30,S31,S32,S33,S40,S41,S42,0,0
int dso::SolidEarthTide::operator()(/*dso::datetime<dso::nanoseconds> &t_tt,
                               double ut1_mjd,*/
                               const Eigen::Matrix<double, 3, 1> &rmoon,
                               const Eigen::Matrix<double, 3, 1> &rsun,
                               std::array<double, 12> &dC,
                               std::array<double, 12> &dS) noexcept {
  // Spherical cordinates of Moon and Sun (ITRF)
  const Eigen::Matrix<double,3,1> Mrfl = dso::car2sph(rmoon);
  const Eigen::Matrix<double,3,1> Srfl = dso::car2sph(rsun);

  // compute associated Legendre functions
  compute_legendre(Mrfl(1), PM.P);
  compute_legendre(Srfl(1), PS.P);

  // Step 1 corrections
  solid_earth_tide_step1(Mrfl(0), Srfl(0), Mrfl(2), Srfl(2), dC, dS);

  // Step 2 corrections
  // double dc20,dc21,ds21,dc22,ds22;
  // solid_earth_tide_step2(dc20,dc21,ds21,dc22,ds22);

  // apply Step 2 corrections
  // dC[0] += dc20;
  // dC[1] += dc21;
  // dS[1] += ds21;
  // dC[2] += dc22;
  // dS[2] += ds22;

  return 0;
}

int dso::SolidEarthTide::acceleration(const Eigen::Matrix<double, 3, 1> &rsat,
                                 const Eigen::Matrix<double, 3, 1> &rmoon,
                                 const Eigen::Matrix<double, 3, 1> &rsun,
                                 Eigen::Matrix<double, 3, 1> &acc) noexcept {
  Eigen::Matrix<double, 3, 1> a;
  Eigen::Matrix<double, 3, 3> partials;

  // compute SH coefficients (C and S)
  std::array<double, 12> dC, dS;
  this->operator()(rmoon, rsun, dC, dS);
  array2harmonics(dC, dS, cs);

  // compute acceleration gradient (ITRF, cartesian)
  test::gravacc3(cs, rsat, degree, cs._Re, cs._GM, a, partials, &V, &W);

  acc = a;

  return 0;
}
