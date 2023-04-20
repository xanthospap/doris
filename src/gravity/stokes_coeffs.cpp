#include "stokes_coeffs.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>

double *dso::StokesCoeffs::allocate(int degree, int order) noexcept {
  deallocate();
  if (degree)
    m_data = new double[(degree + 1) * (degree + 1)];
  m_degree = degree;
  m_order = order;
  return m_data;
}

void dso::StokesCoeffs::deallocate() noexcept {
  if (m_data)
    delete[] m_data;
  m_data = nullptr;
  m_degree = 0;
  m_order = 0;
  return;
}

void dso::StokesCoeffs::resize(int degree, int order) noexcept {
  deallocate();
  allocate(degree, order);
  return;
}

dso::StokesCoeffs::StokesCoeffs(dso::StokesCoeffs &&h) noexcept {
  this->deallocate();
  this->m_degree = h.m_degree;
  this->m_order = h.m_order;
  this->m_data = h.m_data;
  this->_GM = h._GM;
  this->_Re = h._Re;
  this->_cnormalized = h._cnormalized;
  h.m_degree = 0;
  h.m_order = 0;
  h.m_data = nullptr;
}

dso::StokesCoeffs &
dso::StokesCoeffs::operator=(dso::StokesCoeffs &&h) noexcept {
  this->m_degree = h.m_degree;
  this->m_order = h.m_order;
  this->m_data = h.m_data;
  this->_GM = h._GM;
  this->_Re = h._Re;
  this->_cnormalized = h._cnormalized;
  h.m_degree = 0;
  h.m_order = 0;
  h.m_data = nullptr;
  return *this;
}
