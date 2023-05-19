#include "stokes_coeffs.hpp"
#include <stdexcept>

void dso::StokesCoeffs::resize(int degree, int order) noexcept {
  Cnm.resize(degree+1, order+1);
  Snm.resize(degree+1, order+1);
  m_degree = degree;
  m_order = order;
}
  
dso::StokesCoeffs& dso::StokesCoeffs::operator+=(const dso::StokesCoeffs &sc)
{
  if (m_degree<sc.max_degree() || m_order<sc.max_order()) {
    throw std::runtime_error("[ERROR] Failed to add Stokes coefficients\n");
  }

  for (int m=0; m<=sc.max_order(); m++) {
    for (int n=m; n<=sc.max_degree(); n++) {
      C(n,m) += sc.C(n,m);
    }
  }
  
  for (int m=0; m<=sc.max_order(); m++) {
    for (int n=m; n<=sc.max_degree(); n++) {
      S(n,m) += sc.S(n,m);
    }
  }

  return *this;
}
