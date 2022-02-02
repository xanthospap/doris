#ifndef __EARTH_GRAVITY_N_POTENTIAL_HPP__
#define __EARTH_GRAVITY_N_POTENTIAL_HPP__

#include "compact2dmat.hpp"
#include "harmonic_coeffs.hpp"

#ifdef DEBUG
int lagrange_polynomials(double x, double y, double z, double R, int l, int k,
                         Mat2D<MatrixStorageType::RowWise> &V,
                         Mat2D<MatrixStorageType::RowWise> &W) noexcept;
#endif
int lagrange_polynomials(double x, double y, double z, double R, int l, int k,
                         Mat2D<MatrixStorageType::Trapezoid> &V,
                         Mat2D<MatrixStorageType::Trapezoid> &W) noexcept;

#endif