#ifndef __STOKES_HARMONIC_POTENTIAL_COEFFICIENTS_HPP__
#define __STOKES_HARMONIC_POTENTIAL_COEFFICIENTS_HPP__

#ifdef DORIS_EXTRA_CHECKS
#include <cassert>
#include <cstdio>
#endif
#include "cmat2d.hpp"
#include <algorithm>
#include <cstring>

namespace dso {

class StokesCoeffs {
private:
  /* gravitational constant times mass [kg^3/m^2] */
  double _GM{0e0}; 
  /* reference radius of the spherical harmonics [m] */
  double _Re{0e0}; 
  /* coefficients are normaliized (?) */
  bool _cnormalized{true}; 
  /* maximum degree */
  int m_degree{0};
  /* maximum order */
  int m_order{0};
  /* Cnm coefficients */
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> Cnm;
  /* Snm coefficients */
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> Snm;

public:
  StokesCoeffs() noexcept
      : _GM(0e0), _Re(0e0), _cnormalized(true), m_degree(0),
        m_order(0), Cnm{0, 0}, Snm{0, 0} {}

  StokesCoeffs(int n, int m, double GM, double Re)
      : _GM(GM), _Re(Re), _cnormalized(true), m_degree(n), m_order(m),
        Cnm(n + 1, m + 1), Snm(n + 1, m + 1) {}

  StokesCoeffs(int n)
      : _GM(0e0), _Re(0e0), _cnormalized(true), m_degree(n), m_order(n),
        Cnm(n + 1, n + 1), Snm(n + 1, n + 1) {}

  /*
  StokesCoeffs(const StokesCoeffs &h) = delete;
  StokesCoeffs &operator=(const StokesCoeffs &h) = delete;
  StokesCoeffs(StokesCoeffs &&h) noexcept;
  StokesCoeffs &operator=(StokesCoeffs &&h) noexcept;
  ~StokesCoeffs() noexcept { deallocate(); }
  */

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
  double J2() const noexcept { return -Cnm(2,0); };
  void clear() noexcept {
    if (m_degree) {
      Cnm.fill_with(0e0);
      Snm.fill_with(0e0);
    }
  }
  void scale(double factor) noexcept {
    if (m_degree) {
    Cnm.multiply(factor);
    Snm.multiply(factor);
    }
  }
  double C(int n, int m) const noexcept {return Cnm(n,m);}
  double &C(int n, int m) noexcept {return Cnm(n,m);}
  double S(int n, int m) const noexcept {return Snm(n,m);}
  double &S(int n, int m) noexcept {return Snm(n,m);}
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &Cmat() noexcept {return Cnm;}
  dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> &Smat() noexcept {return Snm;}

}; // StokesCoeffs

} // namespace dso

#endif
