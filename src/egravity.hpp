#ifndef __EARTH_GRAVITY_N_POTENTIAL_HPP__
#define __EARTH_GRAVITY_N_POTENTIAL_HPP__

#include "cmat2d.hpp"
#include "eigen3/Eigen/Eigen"
#include "harmonic_coeffs.hpp"
#include "matvec/matvec.hpp"

namespace dso {

/// @brief Computes the perturbational acceleration due to a point mass
/// E.g. use this function we can compute the perturbing acceleration affecting
/// a satellite via sun or moon.
/// @param[in] rsat Satellite position vector
/// @param[in] robj Point mass position vector (e.g. moon)
/// @return A vector containing the acceleration components
/// @see e.g. Curtis, Chapter 10.10
inline Vector3 point_mass_accel(const Vector3 &rsat, const Vector3 &robj,
                                double GM) noexcept {
  //  Relative position vector of satellite w.r.t. point mass
  auto d = rsat - robj;
  // Acceleration
  return -GM * (d / std::pow(d.norm(), 3) + robj / std::pow(robj.norm(), 3));
}

Eigen::Matrix<double, 3, 1>
point_mass_accel(double GM, const Eigen::Matrix<double, 3, 1> &rsat,
                 const Eigen::Matrix<double, 3, 1> &robj,
                 Eigen::Matrix<double, 3, 3> &partials) noexcept;

#ifdef DEBUG
int lagrange_polynomials(
    double x, double y, double z, double R, int l, int k,
    dso::Mat2D<dso::MatrixStorageType::RowWise> &V,
    dso::Mat2D<dso::MatrixStorageType::RowWise> &W) noexcept;
#endif

/// Compute Lagrange polynomials (for spherical harmonics) given a (cartesian)
/// position vector.
/// @param[in] x X-component of position vector in meters [m]
/// @param[in] y Y-component of position vector in meters [m]
/// @param[in] z Z-component of position vector in meters [m]
/// @param[in] Re Earth radius (depending on gravity model)
/// @param[in] l max degree
/// @param[in] k max order (k <= l)
/// @param[out] V Computed values of V lagrange polynomials
/// @param[out] W Computed values for W lagrange polynomials
///
/// Note that for Spherical Harmonics up to C_nn and S_nn, we need to
/// compute the V and W values up to n+1. If we want partials, we must expand
/// to n+2
/// @ref Montenbruck, Gill, Satellite Orbits, Models Methods Applications;
///      ch. 3.2.4, p. 66
int lagrange_polynomials(
    double x, double y, double z, double Re, int l, int k,
    dso::Mat2D<dso::MatrixStorageType::Trapezoid> &V,
    dso::Mat2D<dso::MatrixStorageType::Trapezoid> &W) noexcept;

inline int lagrange_polynomials(
    const Eigen::Matrix<double, 3, 1> &xyz, double Re, int l, int k,
    dso::Mat2D<dso::MatrixStorageType::Trapezoid> &V,
    dso::Mat2D<dso::MatrixStorageType::Trapezoid> &W) noexcept {
  const double x = xyz(0);
  const double y = xyz(1);
  const double z = xyz(2);
  return lagrange_polynomials(x, y, z, Re, l, k, V, W);
}

/// @brief Computes the acceleration due to the harmonic gravity field of the
/// central body
/// @param[in] GM Gravitational coefficient (corresponding to given harmonics)
/// @param[in] Re Reference radius (corresponding to given harmonics)
/// @param[in] hc Spherical harmonics coefficients (un-normalized)
/// @param[in] degree Maximum degree; less or equal to the degree of the hc
/// @param[in] order Maximum order (m_max<=n_max; m_max=0 for zonals, only)
/// @param[out] V Computed values of V lagrange polynomials (computed in an
///               Earth-fixed coordinate system)
/// @param[out] W Computed values for W lagrange polynomials (computed in an
///               Earth-fixed coordinate system)
/// @param[out] acc Acceleration (a=d^2r/dt^2) in x, y, z components in an
///               Earth-fixed coordinate system (the same as the one used
///               to compute the lagrange polynomial values stored in V and W)
/// @warning Note that if we need to compute the potential, aka the harmonic
///          expansion of degree N and order M, then we need the V and W values
///          for degree N+1 and order M+1. Hence, when allocating the structs,
///          users should do something like:
///          Mat2D<MatrixStorageType::Trapezoid> V(degree + 2, order + 2);
///          See Montenbruck 3.2.5 and the example program test_gravacc.cpp.
/// @ref Montenbruck, Gill, Satellite Orbits, Models Methods Applications;
///      ch. 3.2.5, p. 68
int grav_potential_accel(int degree, int order, double Re, double GM,
                         const dso::Mat2D<MatrixStorageType::Trapezoid> &V,
                         const dso::Mat2D<MatrixStorageType::Trapezoid> &W,
                         const dso::HarmonicCoeffs &hc, double *acc) noexcept;

// Note that if we want the partials, V and W indexes must span [0,degree+2]
// (hence the actual number should be degree+3)
int grav_potential_accel(int degree, int order, double Re, double GM,
                         const dso::Mat2D<MatrixStorageType::Trapezoid> &V,
                         const dso::Mat2D<MatrixStorageType::Trapezoid> &W,
                         const dso::HarmonicCoeffs &hc, double *acc,
                         Eigen::Matrix<double, 3, 3> &partials) noexcept;

// Here, the values of Re and GM are extracted from the Gravity model, aka
// the passed in hc instance
inline int
grav_potential_accel(int degree, int order,
                     const dso::Mat2D<MatrixStorageType::Trapezoid> &V,
                     const dso::Mat2D<MatrixStorageType::Trapezoid> &W,
                     const dso::HarmonicCoeffs &hc, double *acc) noexcept {
  return grav_potential_accel(degree, order, hc.Re(), hc.GM(), V, W, hc, acc);
}

// Here, the values of Re and GM are extracted from the Gravity model, aka
// the passed in hc instance.
// Note that if we want the partials, V and W indexes must span [0,degree+2]
// (hence the actual number should be degree+3).
inline int
grav_potential_accel(int degree, int order,
                     const dso::Mat2D<MatrixStorageType::Trapezoid> &V,
                     const dso::Mat2D<MatrixStorageType::Trapezoid> &W,
                     const dso::HarmonicCoeffs &hc, double *acc,
                     Eigen::Matrix<double, 3, 3> &partials) noexcept {
  return grav_potential_accel(degree, order, hc.Re(), hc.GM(), V, W, hc, acc,
                              partials);
}
} // namespace dso

#endif
