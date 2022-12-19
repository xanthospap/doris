#include "associated_legendre.hpp"
#include <cmath>
#include "factorial.hpp"

/// @brief Compute coefficients Anm and Bnm for computation of (normalized)
/// associated Legendre function, using the standard forward column method.
/// When called, it will fill in the *F1 and *F2 coefficient matrices for the 
/// calling instance. Note that this function needs only to be called once, 
/// at instance initialization; it is independent of any spatial/temporal 
/// input.
///
/// @note To be used with the standard forward column method(s). The 
/// Condon–Shortley phase factor (aka (-1)^m) is NOT part of the computation
///
/// @see Holmes et al, 2002
void dso::AssociatedLegendreFunctions::compute_factors() noexcept {
  const int maxDegree = m_degree;
  
  // diagonal terms (Eq. 13)
  F1->operator()(1,1) = std::sqrt(3e0);
  for (int n=2; n<=maxDegree; n++) {
    (*F1)(n,n) = std::sqrt( (2*n+1) / (2e0*n) );
  }

  // compute Anm and Bnm, (Eq. 12)
  for (int m=0; m<=maxDegree; m++) {
    for (int n=m+1; n<=maxDegree; n++) {
      const double f = (2 * n + 1) / (double)((n + m) * (n - m));
      F1->operator()(n,m) = std::sqrt(f*(2*n-1));
      F2->operator()(n,m) = std::sqrt(f*((n-m-1)*(n+m-1))/(2*n-3));
    }
  }

  return;
}

/// @brief Compute the Normalized Associative Legendre Functions, using the 
/// standard forward column method (see Holmes et al, 2002).
/// The function will compute the values of the Normalized P_nm values and 
/// store them in the P coefficient matrix (of the instance).
/// Note that the Condon–Shortley phase factor (aka (-1)^m) is NOT part of the 
/// computation
/// @warning We are assuming here that the coefficient matrices F1 and F2 
/// are already computed.
/// @param[in] angle An angle in [rad], usually the polar angle, colatitude or
///                  geocentric latitude, θ
/// @see Holmes et al, 2002
void dso::AssociatedLegendreFunctions::compute(double angle) noexcept {
  const int maxDegree = m_degree;
  const double u = std::sin(angle);
  const double t = std::cos(angle);

  // for the usage of this scale factor, see Holmes, Sec. 2.7
  const double scaleFactor = 1e280;
  const double scaleFactorInv = 1e-280;

  P(0,0) = 1e0 * scaleFactor;
  
  { // fill elements for order = 0 (Eq. 11)
    const int m = 0;
    P(1,0) = (*F1)(1,0) * P(0,0) * t;
    for (int n=2; n<=maxDegree; n++)
      P(n,m) = (*F1)(n,m) * t * P(n-1, m) - (*F2)(n,m) * P(n-2, m);
  }

  // column wise visiting, for orders: 1<=m<maxDegree, Eq. 11 and 13
  double lastDiagonal = P(0,0);
  for (int m=1; m<maxDegree; m++) {
    // diagonal elements
    P(m,m) = (*F1)(m,m) * u * /*P(m-1,m-1)*/lastDiagonal;
    lastDiagonal = P(m,m);
    // sub-diagonal points
    P(m+1,m) = (*F1)(m+1,m) * t * P(m,m);
    // all other points of same order = m
    for (int n=m+2; n<=maxDegree; n++)
      P(n,m) = (*F1)(n,m) * t * P(n-1, m) - (*F2)(n,m) * P(n-2, m);
  }

  { // we are missing the last (diagonal) element!
    const int m = maxDegree;
    P(m, m) = (*F1)(m,m) * u * /*P(m-1,m-1)*/lastDiagonal;
  }

  // scale
  P.multiply(scaleFactorInv);

  // all done!
  return;
}

/// @brief Compute the Normalized Associative Legendre Functions, using 
/// a semi-analytical approach, see 
/// https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Reparameterization_in_terms_of_angles
/// The function will compute the values of the Normalized P_nm values and 
/// store them in the P coefficient matrix (of the instance).
/// Note that the Condon–Shortley phase factor (aka (-1)^m) is NOT part of the 
/// computation
/// @param[in] angle An angle in [rad], usually the polar angle, colatitude or
///                  geocentric latitude, θ
void dso::AssociatedLegendreFunctions::assign(double angle) noexcept {
  const double u = std::sin(angle);
  const double t = std::cos(angle);
  const double u2 = u * u;
  const double t2 = t * t;

  P(0, 0) = 1e0;
  P(1, 0) = t * std::sqrt(3e0);
  P(2, 0) = .5e0 * std::sqrt(5e0) * (3e0 * t2 - 1e0);
  P(3, 0) = (.5e0 * t * std::sqrt(7e0)) * (5e0 * t2 - 3e0);
  P(4, 0) = (3e0 / 8e0) * t2 * (-30e0 + 35e0 * t2) + (3e0 / 8e0) * 3e0;

  P(1, 1) = u * std::sqrt(3);
  P(2, 1) = 3e0 * std::sqrt(5e0 / 3e0) * t * u;
  P(3, 1) = (3e0 / 2e0) * (5e0 * t2 - 1e0) * u * std::sqrt((7e0) / 6e0);
  P(4, 1) = (2.5e0*std::sqrt(9e0/10e0)) * u * ((7e0 * t2 - 3e0) * t);


  P(2, 2) = 3e0 * std::sqrt(5e0 / 12e0) * u2;
  P(3, 2) = 15e0 * std::sqrt(7e0 / 60e0) * t * u2;
  P(4, 2) = (7.5e0 * std::sqrt(5e-2)) * (7e0 * t2 - 1e0) * u2;

  P(3, 3) = 15e0 * u2 * u * std::sqrt(14e0 / 720e0);
  P(4, 3) = 105e0 * std::sqrt(18e0 / 5040e0) * t * u2 * u;

  P(4, 4) = std::sqrt(18e0 / 40320e0) * 105e0 * u2 * u2;

  return;
}
