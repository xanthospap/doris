#ifndef __HARMONIC_POTENTIAL_COEFFICIENTS_HPP__
#define __HARMONIC_POTENTIAL_COEFFICIENTS_HPP__

/// Define a data structure to hold Harmonics coefficients. This structure
/// stores a two-dimensional data matrix, where each row corresponds to the
/// S and C coefficients of a certain degree. To save space, the data matrix
/// is in compact form, as follows (degree=n, order=m):
/// C_00 S_n1 S_n2     S_n3     ... S_nm
/// C_10 C_11 S_(n-1)1 S_(n-1)2 ... S_(n-1)m
/// ...
/// C_(n-1)0 C_(n-1)1 C_(n-1)2 .... C_(n-1)m S_11
/// C_n0     C_n1     C_n2     .... C_n(m-1) C_nm
/// If degree > order, empty space is left between the S and C coefficients.

#ifdef DEBUG
#include <cassert>
#include <cstdio>
#endif

namespace dso {

/// @todo Why the fuck did i make this a pointer to pointer? Could be just
/// one big chunk of memory.
///
/// @brief Storage and access of Harmonic Coefficients.
/// We are allocating and using a 2-d array of size (degree+1) * (degree+1) 
/// (hence not actually using the order of the Coefficients). This make take 
/// up more space than actually needed but is a bit faster.
/// The struct is 'agnostic' concerning the order of the coefficients. It is
/// the user's responsibility to use consistent values for order.
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
class HarmonicCoeffs {
public:
  double _GM{0e0}; ///< gravitational constant times mass of the earth
  double _Re{0e0}; ///< reference radius of the spherical harmonic development
  bool _cnormalized{true};  ///< coefficients are normaliized (?)
  int m_degree{0};          ///< maximum degree
  double **m_data{nullptr}; ///< the actual data/coefficients

  /// @brief allocate memory to hold the data.
  double *allocate() noexcept;

  /// @brief free memmory used by the structure.
  int deallocate() noexcept;

public:
  HarmonicCoeffs()
      : _GM(0e0), _Re(0e0), _cnormalized(true), m_degree(0), m_data(nullptr){};

  HarmonicCoeffs(int n, double GM, double Re)
      : _GM(GM), _Re(Re), _cnormalized(true), m_degree(n) {
    allocate();
  }

  HarmonicCoeffs(int n) : _GM(0e0), _Re(0e0), _cnormalized(true), m_degree(n) {
    allocate();
  }

  HarmonicCoeffs(const HarmonicCoeffs &h) = delete;

  HarmonicCoeffs &operator=(const HarmonicCoeffs &h) = delete;

  HarmonicCoeffs(HarmonicCoeffs &&h) noexcept;

  HarmonicCoeffs &operator=(HarmonicCoeffs &&h) noexcept;

  ~HarmonicCoeffs() noexcept { deallocate(); }

  /// @brief Resize; check current capacity and only re-allocated data if 
  ///        needed. m_degree set to new value.
  void resize(int degree) noexcept;

#ifdef DEBUG
  void print(double scale = 1e0) noexcept {
    for (int i = 0; i <= m_degree; i++) {
      for (int j = 0; j <= m_degree; j++) {
        printf("%15.10e ", scale * m_data[i][j]);
      }
      printf("\n");
    }
  }
#endif

  int degree() const noexcept { return m_degree; }
  double GM() const noexcept { return _GM; }
  double Re() const noexcept { return _Re; }
  bool normalized() const noexcept { return _cnormalized; }
  double &GM() noexcept { return _GM; }
  double &Re() noexcept { return _Re; }
  bool &normalized() noexcept { return _cnormalized; }
  double J2() const noexcept {
#ifdef DEBUG
    assert( m_data[2][0]  == this->C(2,0));
#endif
    return -m_data[2][0];
  };

  /// @brief De-normalize harmonic coefficients.
  /// If the coefficients are alredy un-normalized, this is a no-op.
  /// @ref Montenbruck, Gill, Satellite Orbits, Models Methods Applications;
  /// See Eq. 3.13 in Chapter 3.2
  int denormalize(int order = -1) noexcept;

  /// @brief Get a pointer to the C coefficients of degree 'degree'.
  /// E.g. C_row[5] will hold the C(5,0) coefficient, C_row[5] + 1 will point
  /// to the C(5,1) coefficient and C_row[5] + 5 will point to the C(5,5)
  /// coefficient (that is degree=5 and order=5).
  double *C_row(int degree) noexcept {
#ifdef DEBUG
    assert(degree <= m_degree);
#endif
    return m_data[degree]; // C(degree,0)-> C(degree, degree)
  }

  /// @brief Get a pointer to the C coefficients of degree 'degree'.
  /// E.g. C_row[5] will hold the C(5,0) coefficient, C_row[5] + 1 will point
  /// to the C(5,1) coefficient and C_row[5] + 5 will point to the C(5,5)
  /// coefficient (that is degree=5 and order=5).
  const double *C_row(int degree) const noexcept {
#ifdef DEBUG
    assert(degree <= m_degree);
#endif
    return m_data[degree]; // C(degree,0)-> C(degree, degree)
  }

  /// @brief Get the C coefficient of degree i and order j
  double &C(int i, int j) noexcept { 
#ifdef DEBUG
    assert(i <= m_degree && j <= i);
#endif
    return C_row(i)[j];
  }

  /// @brief Get the C coefficient of degree i and order j
  const double &C(int i, int j) const noexcept { 
#ifdef DEBUG
    assert(i <= m_degree && j <= i);
#endif
    return C_row(i)[j];
  }

  /// @brief Get a pointer to the S coefficients of degree 'degree'.
  /// @warning Sn0 (aka S coefficients for degree=0) are always zero and are
  ///          not stored in the data array. Hence, in contrast to the
  ///          corresponding S_row() function, this function will not return a
  ///          pointer to S_0m, a pointer to S_1m.
  /// E.g. S_row[5] will hold the S(5,1) coefficient, C_row[5] + 1 will point
  /// to the C(5,2) coefficient and S_row[5] + 5 will point to the S(5,5)
  /// coefficient (that is degree=5 and order=5).
  double *S_row(int degree) noexcept {
#ifdef DEBUG
    assert(degree <= m_degree && degree != 0);
#endif
    int off = m_degree - degree;
    return m_data[off] + off + 1; // S(degree,1)-> S(degree, degree)
  }

  /// @brief Get a pointer to the S coefficients of degree 'degree'.
  /// @warning Sn0 (aka S coefficients for degree=0) are always zero and are
  ///          not stored in the data array. Hence, in contrast to the
  ///          corresponding S_row() function, this function will not return a
  ///          pointer to S_0m, a pointer to S_1m.
  /// E.g. S_row[5] will hold the S(5,1) coefficient, C_row[5] + 1 will point to
  /// the C(5,2) coefficient and S_row[5] + 5 will point to the S(5,5)
  /// coefficient (that is degree=5 and order=5).
  const double *S_row(int degree) const noexcept {
#ifdef DEBUG
    assert(degree <= m_degree && degree != 0);
#endif
    int off = m_degree - degree;
    return m_data[off] + off + 1; // S(degree,1)-> S(degree, degree)
  }

  /// @brief Get the S coefficient of degree i and order j
  /// @warning never ask for a coefficient with order (aka j) = 0. All S_nm for
  ///          m=0 are equal to 0e0
  double &S(int i, int j) noexcept { 
#ifdef DEBUG
    assert(i <= m_degree && j <= i && i!=0);
#endif
    return S_row(i)[j - 1];
  }

  /// @brief Get the S coefficient of degree i and order j
  /// @warning never ask for a coefficient with order (aka j) = 0. All S_nm for
  ///          m=0 are equal to 0e0
  const double &S(int i, int j) const noexcept { 
#ifdef DEBUG
    assert(i <= m_degree && j <= i && i!=0);
#endif
    return S_row(i)[j - 1];
  }

}; // HarmonicCoeffs

} // dso

#endif
