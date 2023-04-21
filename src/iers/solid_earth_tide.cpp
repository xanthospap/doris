#include "egravity.hpp"
#include "geodesy/geodesy.hpp"
#include "stokes_coeffs.hpp"
#include "iers2010/iau.hpp"
#include "iers2010/iersc.hpp"
#include "tides.hpp"
#ifdef DEBUG
#include "geodesy/units.hpp" // only for debugging
#endif

namespace {
/// @param[in] dC   Normalized corrections to C coefficients, in the order:
///            dC = C20,C21,C22,C30,C31,C32,C33,C40,C41,C42,0,0
///                   0   1   2   3   4   5   6   7   8   9
/// @param[in] dS   Normalized corrections to S coefficients, in the order:
///            dS = 0,S21,S22,0,S31,S32,S33,0,S41,S42,0,0
///                 0,  1   2 3   4   5   6 7   8   9
int array2harmonics(const std::array<double, 12> &dC,
                    const std::array<double, 12> &dS,
                    dso::StokesCoeffs &cs) noexcept {
  cs.clear();
  cs.C(2, 0) = dC[0];
  cs.C(3, 0) = dC[3];
  cs.C(4, 0) = dC[7];
  cs.C(2, 1) = dC[1];
  cs.C(3, 1) = dC[4];
  cs.C(4, 1) = dC[8];
  cs.C(2, 2) = dC[2];
  cs.C(3, 2) = dC[5];
  cs.C(4, 2) = dC[9];
  cs.C(3, 3) = dC[6];

  cs.S(2, 1) = dS[1];
  cs.S(2, 2) = dS[2];
  cs.S(3, 1) = dS[4];
  cs.S(3, 2) = dS[5];
  cs.S(3, 3) = dS[6];
  cs.S(4, 1) = dS[8];
  cs.S(4, 2) = dS[9];

  return 0;
}
} // unnamed namespace

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
int dso::SolidEarthTide::operator()(const dso::TwoPartDate &mjdtt,
                                    [[maybe_unused]]const dso::TwoPartDate &mjdut1,
                                    const Eigen::Matrix<double, 3, 1> &rmoon,
                                    const Eigen::Matrix<double, 3, 1> &rsun,
                                    std::array<double, 12> &dC,
                                    std::array<double, 12> &dS) noexcept {
  
  // Step 1 corrections, get δC and δC
  solid_earth_tide_step1(rmoon, rsun, dC, dS);

  // Step 2 corrections
  double dc20, dc21, ds21, dc22, ds22;
  solid_earth_tide_step2(mjdtt, mjdut1, dc20, dc21, ds21, dc22, ds22);

  // apply Step 2 corrections (to Step 1)
  dC[0] += dc20;
  dC[1] += dc21;
  dS[1] += ds21;
  dC[2] += dc22;
  dS[2] += ds22;

  return 0;
}

int dso::SolidEarthTide::acceleration(
    const dso::TwoPartDate &mjdtt, const dso::TwoPartDate &mjdut1,
    const Eigen::Matrix<double, 3, 1> &rsat,
    const Eigen::Matrix<double, 3, 1> &rmoon,
    const Eigen::Matrix<double, 3, 1> &rsun, Eigen::Matrix<double, 3, 1> &acc,
    Eigen::Matrix<double, 3, 3> &acc_gradient) noexcept {

  Eigen::Matrix<double, 3, 1> a = Eigen::Matrix<double, 3, 1>::Zero();
  acc_gradient = Eigen::Matrix<double,3,3>::Zero();

  // compute SH coefficient corrections (δC and δS) for Step1 and Step2
  std::array<double, 12> dC, dS;
  this->operator()(mjdtt, mjdut1, rmoon, rsun, dC, dS);

  // transform array to harmonics
  array2harmonics(dC, dS, cs);

  // compute acceleration at satellite position (ITRF, cartesian)
  dso::gravity_acceleration(cs, rsat, degree, cs.Re(), cs.GM(), a, acc_gradient, &V, &W);

  acc = a;

  return 0;
}

/*
int dso::SolidEarthTide::acceleration(
    const dso::TwoPartDate &mjdtt, const dso::TwoPartDate &mjdut1,
    double xp_sec, double yp_sec,
    const Eigen::Matrix<double, 3, 1> &rsat,
    const Eigen::Matrix<double, 3, 1> &rmoon,
    const Eigen::Matrix<double, 3, 1> &rsun, Eigen::Matrix<double, 3, 1> &acc,
    Eigen::Matrix<double, 3, 3> &acc_gradient) noexcept {

  Eigen::Matrix<double, 3, 1> a = Eigen::Matrix<double, 3, 1>::Zero();
  acc_gradient = Eigen::Matrix<double,3,3>::Zero();

  // compute SH coefficient corrections (δC and δS) for Step1 and Step2
  std::array<double, 12> dC, dS;
  this->operator()(mjdtt, mjdut1, rmoon, rsun, dC, dS);

  // compute Pole Tide
  //auto [dc21,ds21] = SolidEarthPoleTide().poleTide(mjdtt, xp_sec, yp_sec);
  //dC[1] += dc21;
  //dS[1] += ds21;

  // transform array to harmonics
  array2harmonics(dC, dS, cs);

  // compute acceleration at satellite position (ITRF, cartesian)
  dso::gravity_acceleration(cs, rsat, degree, cs.Re(), cs.GM(), a, acc_gradient, &V, &W);

  acc = a;

  return 0;
}
*/
