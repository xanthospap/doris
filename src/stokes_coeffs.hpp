#ifndef __STOKES_HARMONIC_POTENTIAL_COEFFICIENTS_HPP__
#define __STOKES_HARMONIC_POTENTIAL_COEFFICIENTS_HPP__

#ifdef DORIS_EXTRA_CHECKS
#include <cassert>
#include <cstdio>
#endif
#include <algorithm>
#include <cstring>

namespace dso {

/// As mentioned, internaly the sctructure holds coefficients spanning S/C
/// from [0-n] and [0-m] for degree and order respectively. Hence, the
/// following is perfectly legal and what **should** be done:
/// int degree = 20;
/// int order  = 15;
/// HarmonicCoeffs HC(degree);
/// ....
/// HC.C(20,15) = ....
/// Aka, you can ask for the coefficient with degree=20 and order=15.
///
/// However, you cannot ask for S coefficients with order=0. These are always
/// equal to zero and **are not stored** in the structure.
///
/// @note No checks are performed for the validity of the degree/order indexes
///       when asking the structure to privide an element (or a row of C or S
///       coefficients. It is the user's responsibility to use correct values.
///
/// @warning Snm coefficients are always zero for m=0. Hence, these values are
///       not even stored in the storage matrix and asking for them, e.g.
///       HarmonicCoeffs hc(...);
///       ....
///       hc.S(1,0); // Error!!
class StokesCoeffs {
private:
  double _GM{0e0}; /* gravitational constant times mass [kg^3/m^2] */
  double _Re{0e0}; /* reference radius of the spherical harmonics [m] */
  bool _cnormalized{true}; /* coefficients are normaliized (?) */
  int m_degree{0};         /* maximum degree */
  int m_order{0};          /* maximum order */
  double *m_data{nullptr}; /* the actual data/coefficients */

  /* @brief allocate memory to hold the data */
  double *allocate(int degree, int order) noexcept;

  /* @brief free memmory used by the structure */
  void deallocate() noexcept;

public:
  StokesCoeffs() noexcept
      : _GM(0e0), _Re(0e0), _cnormalized(true), m_degree(0), m_order(0),
        m_data(nullptr){};

  StokesCoeffs(int n, int m, double GM, double Re)
      : _GM(GM), _Re(Re), _cnormalized(true), m_degree(n), m_order(m) {
    allocate(n, m);
  }

  StokesCoeffs(int n)
      : _GM(0e0), _Re(0e0), _cnormalized(true), m_degree(n), m_order(n) {
    allocate(n, n);
  }

  StokesCoeffs(const StokesCoeffs &h) = delete;
  StokesCoeffs &operator=(const StokesCoeffs &h) = delete;
  StokesCoeffs(StokesCoeffs &&h) noexcept;
  StokesCoeffs &operator=(StokesCoeffs &&h) noexcept;
  ~StokesCoeffs() noexcept { deallocate(); }

  /* @brief Resize; check current capacity and only re-allocated data if
   *      needed. m_degree set to new value.
   */
  void resize(int degree, int order) noexcept;

  int max_degree() const noexcept { return m_degree; }
  int max_order() const noexcept { return m_order; }
  double GM() const noexcept { return _GM; }
  double Re() const noexcept { return _Re; }
  bool normalized() const noexcept { return _cnormalized; }
  double &GM() noexcept { return _GM; }
  double &Re() noexcept { return _Re; }
  bool &normalized() noexcept { return _cnormalized; }
  double J2() const noexcept { return -m_data[2 * (m_degree + 1)]; };
  void clear() noexcept {
    if (m_data)
      std::memset(m_data, 0, (m_degree + 1) * (m_degree + 1) * sizeof(double));
  }
  void scale(double factor) noexcept {
    std::for_each(m_data, m_data + (m_degree + 1) * (m_degree + 1),
                  [=](double &d) { d *= factor; });
  }

  /// @brief De-normalize harmonic coefficients.
  /// If the coefficients are alredy un-normalized, this is a no-op.
  /// @ref Montenbruck, Gill, Satellite Orbits, Models Methods Applications;
  /// See Eq. 3.13 in Chapter 3.2
  /// TODO this needs to be better, see for example
  /// AssociatedLegendreFunctions::normalize
  int denormalize(int order = -1) noexcept;

  /// @brief Get a pointer to the C coefficients of order 'order'.
  /// E.g. C_row[5] will hold the C(5,0) coefficient, C_row[5] + 1 will point
  /// to the C(5,1) coefficient and C_row[5] + 5 will point to the C(5,5)
  /// coefficient (that is degree=5 and order=5).
  double *C_col(int order) noexcept {
#ifdef DORIS_EXTRA_CHECKS
    assert(order <= m_order);
#endif
    return m_data + (m_degree+1)*order; /* C(order,order) to C(degree,order) */
  }

  /// @brief Get a pointer to the C coefficients of degree 'degree'.
  /// E.g. C_row[5] will hold the C(5,0) coefficient, C_row[5] + 1 will point
  /// to the C(5,1) coefficient and C_row[5] + 5 will point to the C(5,5)
  /// coefficient (that is degree=5 and order=5).
  const double *C_col(int order) const noexcept {
#ifdef DORIS_EXTRA_CHECKS
    assert(order <= m_order);
#endif
    /* C(order,order) to C(degree,order) */
    return m_data + (m_degree+1)*order;
  }

  /// @brief Get the C coefficient of degree i and order j
  double &C(int n, int m) noexcept {
#ifdef DORIS_EXTRA_CHECKS
    assert(n <= m_degree && m <= n);
#endif
    return C_col(m)[n-m];
  }

  /// @brief Get the C coefficient of degree i and order j
  const double &C(int n, int m) const noexcept {
#ifdef DORIS_EXTRA_CHECKS
    assert(n <= m_degree && m <= n);
#endif
    return C_col(m)[n-m];
  }

  /// @brief Get a pointer to the S coefficients of degree 'degree'.
  /// @warning Sn0 (aka S coefficients for degree=0) are always zero and are
  ///          not stored in the data array. Hence, in contrast to the
  ///          corresponding S_row() function, this function will not return a
  ///          pointer to S_0m, a pointer to S_1m.
  /// E.g. S_row[5] will hold the S(5,1) coefficient, C_row[5] + 1 will point
  /// to the C(5,2) coefficient and S_row[5] + 5 will point to the S(5,5)
  /// coefficient (that is degree=5 and order=5).
  double *S_col(int order) noexcept {
#ifdef DORIS_EXTRA_CHECKS
    assert(order <= m_order && order != 0);
#endif
    /* S(order,order) to S(n,order), for order>1 */
    return m_data + (m_degree+1)*(m_degree+1)-order*m_degree;
  }

  /// @brief Get a pointer to the S coefficients of degree 'degree'.
  /// @warning Sn0 (aka S coefficients for degree=0) are always zero and are
  ///          not stored in the data array. Hence, in contrast to the
  ///          corresponding S_row() function, this function will not return a
  ///          pointer to S_0m, a pointer to S_1m.
  /// E.g. S_row[5] will hold the S(5,1) coefficient, C_row[5] + 1 will point to
  /// the C(5,2) coefficient and S_row[5] + 5 will point to the S(5,5)
  /// coefficient (that is degree=5 and order=5).
  const double *S_col(int order) const noexcept {
#ifdef DORIS_EXTRA_CHECKS
    assert(order <= m_order && order != 0);
#endif
    /* S(order,order) to S(n,order), for order>1 */
    return m_data + (m_degree+1)*(m_degree+1)-order*m_degree;
  }

  /// @brief Get the S coefficient of degree i and order j
  /// @warning never ask for a coefficient with order (aka j) = 0. All S_nm for
  ///          m=0 are equal to 0e0
  double &S(int n, int m) noexcept {
#ifdef DORIS_EXTRA_CHECKS
    assert(n <= m_degree && m <= n && m != 0);
#endif
    return S_col(m)[n-m];
  }

  /// @brief Get the S coefficient of degree i and order j
  /// @warning never ask for a coefficient with order (aka j) = 0. All S_nm for
  ///          m=0 are equal to 0e0
  const double &S(int n, int m) const noexcept {
#ifdef DORIS_EXTRA_CHECKS
    assert(n <= m_degree && m <= n && m != 0);
#endif
    return S_col(m)[n-m];
  }

}; // StokesCoeffs

} // namespace dso

#endif
