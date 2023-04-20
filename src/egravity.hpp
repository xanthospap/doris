#ifndef __EARTH_GRAVITY_N_POTENTIAL_HPP__
#define __EARTH_GRAVITY_N_POTENTIAL_HPP__

#include "associated_legendre.hpp"
#include "cmat2d.hpp"
#include "datetime/dtcalendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "stokes_coeffs.hpp"
#include <cassert>
#ifdef DEBUG
#include <cstdio>
#endif

namespace dso {
int gravity_acceleration(const dso::StokesCoeffs &cs,
                         const Eigen::Matrix<double, 3, 1> &p, int degree,
                         double Re, double GM, Eigen::Matrix<double, 3, 1> &acc,
                         Eigen::Matrix<double, 3, 3> &gradient) noexcept;
int gravity_acceleration(
    const dso::StokesCoeffs &cs, const Eigen::Matrix<double, 3, 1> &p,
    int degree, double Re, double GM, Eigen::Matrix<double, 3, 1> &acc,
    Eigen::Matrix<double, 3, 3> &gradient,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *W,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *M) noexcept;

/// @brief Parse a gravity model (given in icgem format) to a HarmonicCoeffs
///        instance
/// @param[in] model_fn An icegm file, describing a static gravity/geopotential
///        model
/// @param[in] degree Max degree for harmonics to extract; should be at
///        maximum equal to the model's max degree
/// @param[in] order Max order for harmonics to extract; should be at
///        maximum equal to the model's max order and <= to degree
/// @param[out] harmonics An dso::HarmonicCoeffs instance to store results. If
///        (at input) the instance has wrong size, it will be resized to match
///        the needs of the requested harmonics
/// @param[in] denormalize Denormalize coefficients (if normalized)
/// @return Anything other than 0 denotes an error
int parse_gravity_model(const char *model_fn, int degree, int order,
                        const dso::TwoPartDate &t,
                        dso::StokesCoeffs &harmonics) noexcept;

/// @brief Computes the perturbational acceleration due to a point mass
/// E.g. use this function we can compute the perturbing acceleration affecting
/// a satellite via sun or moon.
/// @param[in] rsat Satellite position vector
/// @param[in] robj Point mass position vector (e.g. moon)
/// @return A vector containing the acceleration components
/// @see e.g. Curtis, Chapter 10.10
inline Eigen::Matrix<double, 3, 1>
point_mass_accel(const Eigen::Matrix<double, 3, 1> &rsat,
                 const Eigen::Matrix<double, 3, 1> &robj, double GM) noexcept {
  //  Relative position vector of satellite w.r.t. point mass
  auto d = rsat - robj;
  // Acceleration
  return -GM * (d / std::pow(d.norm(), 3) + robj / std::pow(robj.norm(), 3));
}

Eigen::Matrix<double, 3, 1>
point_mass_accel(double GM, const Eigen::Matrix<double, 3, 1> &rsat,
                 const Eigen::Matrix<double, 3, 1> &robj,
                 Eigen::Matrix<double, 3, 3> &partials) noexcept;

Eigen::Matrix<double, 3, 1>
point_mass_accel(double GM, const Eigen::Matrix<double, 3, 1> &rsat,
                 const Eigen::Matrix<double, 3, 1> &robj) noexcept;

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

} // namespace dso

#endif
