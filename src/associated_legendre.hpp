#ifndef __ASSICIATED_LEGENDRE_FUNCTIONS_HPP__
#define __ASSICIATED_LEGENDRE_FUNCTIONS_HPP__

#include "cmat2d.hpp"

namespace dso {
namespace detail {
  constexpr const int MIN_ASSOCIATEDLEGENDREFUNCTIONS_DEGREE = 4;
};//detail
struct AssociatedLegendreFunctions {
  int m_degree; ///< degree, inclusive; order  = degree
  ///< Lower triangular, of size (degree+1, degree+1)
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> P;

  AssociatedLegendreFunctions(int degree)
      : m_degree(degree),
        P((degree < detail::MIN_ASSOCIATEDLEGENDREFUNCTIONS_DEGREE)
              ? (detail::MIN_ASSOCIATEDLEGENDREFUNCTIONS_DEGREE + 1)
              : (degree + 1),
          (degree < detail::MIN_ASSOCIATEDLEGENDREFUNCTIONS_DEGREE)
              ? (detail::MIN_ASSOCIATEDLEGENDREFUNCTIONS_DEGREE + 1)
              : (degree + 1)) {}

  int degree() const noexcept {return m_degree;}

  double &operator()(int n, int m) noexcept {
#ifdef DEBUG
    assert(n>=0 && m <= n && n <= m_degree);
#endif
    return P(n,m);
  }
  
  double operator()(int n, int m) const noexcept {
#ifdef DEBUG
    assert(n>=0 && m <= n && n <= m_degree);
#endif
    return P(n,m);
  }

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
  /// void compute(double phi) noexcept;

  /// @brief Normalize Legendere functions, using the scale factor:
  ///        N_nm = sqrt[ (n-m)! (2n+1) (2-Kronecker(0,m) / (n+m)! ]
  void normalize() noexcept;
}; // AssociatedLegendreFunction
} // namespace dso

#endif
