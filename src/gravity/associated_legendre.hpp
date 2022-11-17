#ifndef __ASSICIATED_LEGENDRE_FUNCTIONS_HPP__
#define __ASSICIATED_LEGENDRE_FUNCTIONS_HPP__

#include "cmat2d.hpp"

namespace dso {
struct AssociatedLegendreFunctions {
  int m_order; ///< order, inclusive
  ///< Lower triangular, of size (order+1, order+1)
  dso::Mat2D<dso::MatrixStorageType::LwTriangularRowWise> P;

  AssociatedLegendreFunctions(int order)
      : m_order(order), P(order + 1, order + 1){};

  /// The computation algorithm follow Montenbruck et al, 2012, Section 3.2.4
  /// Starting with P_00 = 1, all polynomials up to the desired degree and
  /// order are first calculated from:
  /// P_mm(u) = (2m-1) x (1-u**2)^(1/2) x P_(m-1)(m-1) (Eq. 3.23)
  /// where u -> sin(f) and (1-u**2)^(1/2) -> cos(f). With these results,
  /// the remaining values may be obtained from
  /// P_(m+1)(m) = (2m+1) x u x P_mm(u)                (Eq. 3.24)
  /// and from the recursion:
  /// P_nm(u) = (1/(n-m)) [(2n-1) u P_(n-1)(m)(u) - (n+m-1)P_(n-2)(m)(u)]
  /// for n > m + 1
  /// TODO This is more of a column-wise computation algorithm. Could maybe 
  /// gain speed if we used LwTriangularColWise storage type or use another 
  /// algorithm
  void compute(double phi) noexcept;

  /// @brief Normalize Legendere functions, using the scale factor:
  ///        N_nm = sqrt[ (n-m)! (2n+1) (2-Kronecker(0,m) / (n+m)! ]
  void normalize() noexcept;
}; // AssociatedLegendreFunction
} // namespace dso

#endif
