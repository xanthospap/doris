#ifndef __EARTH_GRAVITY_N_POTENTIAL_HPP__
#define __EARTH_GRAVITY_N_POTENTIAL_HPP__

#include "compact2dmat.hpp"
#include "harmonic_coeffs.hpp"
#include "iers2010/matvec.hpp"

namespace dso {

/// @brief Computes the perturbational acceleration due to a point mass
/// E.g. use this function we can compute the perturbing acceleration affecting
/// a satellite via sun or moon.
/// @param[in] rsat Satellite position vector
/// @param[in] robj Point mass position vector (e.g. moon)
/// @return A vector containing the acceleration components
/// @see e.g. Curtis, Chapter 10.10
inline
Vector3 point_mass_accel(const Vector3 &rsat, const Vector3 &robj,
                         double GM) noexcept {
  //  Relative position vector of satellite w.r.t. point mass
  auto d = rsat - robj;
  // Acceleration
  return -GM * (d / std::pow(d.norm(), 3) + robj / std::pow(robj.norm(), 3));
}

#ifdef DEBUG
int lagrange_polynomials(
    double x, double y, double z, double R, int l, int k,
    dso::Mat2D<dso::MatrixStorageType::RowWise> &V,
    dso::Mat2D<dso::MatrixStorageType::RowWise> &W) noexcept;
#endif

/// Compute Lagrange polynomials (for spherical harmonics) given a (cartesian) 
/// position vector.
/// @param[in] x X-component of position vector in meters
/// @param[in] y Y-component of position vector in meters
/// @param[in] z Z-component of position vector in meters
/// @param[in] R Earth radius (depending on gravity model)
/// @param[in] l max degree
/// @param[in] k max order (k <= l)
/// @param[out] V Computed values of V lagrange polynomials
/// @param[out] W Computed values for W lagrange polynomials
/// 
/// @ref Montenbruck, Gill, Satellite Orbits, Models Methods Applications;
///      ch. 3.2.4, p. 66
int lagrange_polynomials(
    double x, double y, double z, double R, int l, int k,
    dso::Mat2D<dso::MatrixStorageType::Trapezoid> &V,
    dso::Mat2D<dso::MatrixStorageType::Trapezoid> &W) noexcept;

/// @brief Computes the acceleration due to the harmonic gravity field of the
/// central body
/// @param[in] GM Gravitational coefficient (corresponding to given harmonics)
/// @param[in] Re Reference radius (corresponding to given harmonics)
/// @param[in] hc Spherical harmonics coefficients (un-normalized)
/// @param[in] degree Maximum degree; less or equal to the degree of the hc
/// @param[in] order Maximum order (m_max<=n_max; m_max=0 for zonals, only)
/// @param[out] V Computed values of V lagrange polynomials
/// @param[out] W Computed values for W lagrange polynomials
/// @param[out] acc Acceleration (a=d^2r/dt^2) in x, y, z components
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
}// dso

#endif
