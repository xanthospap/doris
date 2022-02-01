#include "harmonic_coeffs.hpp"

int HarmonicCoeffs::allocate() noexcept {
    m_data = new double*[m_degree+1];
    for (int i=0; i<=m_degree; i++)
        m_data[i] = new double[m_degree+1];
    return 0;
}

int HarmonicCoeffs::deallocate() noexcept {
  if (m_data) {
    for (int i = 0; i <= m_degree; i++)
      delete[] m_data[i];
    delete [] m_data;
  }
  return 0;
}

HarmonicCoeffs::HarmonicCoeffs(HarmonicCoeffs &&h) noexcept {
  this->deallocate();
  this->m_degree = h.m_degree;
  this->m_data = h.m_data;
  h.m_degree = 0;
  h.m_data = nullptr;
}

HarmonicCoeffs &HarmonicCoeffs::operator=(HarmonicCoeffs &&h) noexcept {
  this->m_degree = h.m_degree;
  this->m_data = h.m_data;
  h.m_degree = 0;
  h.m_data = nullptr;
  return *this;
}