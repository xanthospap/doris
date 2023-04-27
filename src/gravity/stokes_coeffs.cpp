#include "stokes_coeffs.hpp"

void dso::StokesCoeffs::resize(int degree, int order) noexcept {
  Cnm.resize(degree+1, order+1);
  Snm.resize(degree+1, order+1);
  m_degree = degree;
  m_order = order;
}
