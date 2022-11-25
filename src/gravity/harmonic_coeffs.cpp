#include "harmonic_coeffs.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>

double *dso::HarmonicCoeffs::allocate(int degree, int order) noexcept {
  if (m_data)
    delete[] m_data;
  m_data = nullptr;
  if (degree)
    m_data = new double[(degree + 1) * (degree + 1)];
  m_degree = degree;
  m_order = order;
  return m_data;
}

void dso::HarmonicCoeffs::deallocate() noexcept {
  if (m_data)
    delete[] m_data;
  m_data = nullptr;
  m_degree = 0;
  m_order = 0;
  return;
}

void dso::HarmonicCoeffs::resize(int degree, int order) noexcept {
  deallocate();
  allocate(degree, order);
  return;
}

dso::HarmonicCoeffs::HarmonicCoeffs(dso::HarmonicCoeffs &&h) noexcept {
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

dso::HarmonicCoeffs &
dso::HarmonicCoeffs::operator=(dso::HarmonicCoeffs &&h) noexcept {
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

double _sum(int n, int m) noexcept {
  if (m == 0)
    return 1e0;
  int start = n - m + 1;
  int stop = n + m;
  double fac = 1e0;
  for (int i = start; i <= stop; i++)
    fac *= (double)i;
  return fac;
}

int compute_fac(int n, int m, double *fac) noexcept {
  for (int i = 0; i <= m; i++)
    fac[i] = std::sqrt(_sum(n, i) / (2 * n + 1) / (2 - (i == 0)));
  return 0;
}

int dso::HarmonicCoeffs::denormalize(int order) noexcept {
  if (!normalized())
    return 0;

  if (order < 0)
    order = m_degree;
  if (order > m_degree) {
    fprintf(stderr,
            "[ERROR] Invalid order provided (%d) for harmonics of max degree "
            "%d! Cannot denormalize (traceback: %s)\n",
            order, m_degree, __func__);
    return 1;
  }

  assert(m_degree > 0);
  double *facs = new double[m_degree + 1];
  for (int n = 1; n <= m_degree; n++) {
    int mm = std::min(order, n);
    compute_fac(n, mm, facs);

    // C coeffs, 0-std::min(order,n)
    double *c = this->C_row(n);
    for (int m = 0; m <= mm; m++)
      c[m] /= facs[m];

    // S coeffs, 1--std::min(order,n)
    double *s = this->S_row(n);
    for (int m = 1; m <= mm; m++)
      s[m - 1] /= facs[m];
  }
  delete[] facs;

  this->_cnormalized = false;
  return 0;
}
