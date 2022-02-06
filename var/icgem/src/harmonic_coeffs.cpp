#include "harmonic_coeffs.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>

int dso::HarmonicCoeffs::allocate() noexcept {
  m_data = new double *[m_degree + 1];
  for (int i = 0; i <= m_degree; i++)
    m_data[i] = new double[m_degree + 1];
  return 0;
}

int dso::HarmonicCoeffs::deallocate() noexcept {
  if (m_data) {
    for (int i = 0; i <= m_degree; i++)
      delete[] m_data[i];
    delete[] m_data;
  }
  return 0;
}

dso::HarmonicCoeffs::HarmonicCoeffs(dso::HarmonicCoeffs &&h) noexcept {
  this->deallocate();
  this->m_degree = h.m_degree;
  this->m_data = h.m_data;
  h.m_degree = 0;
  h.m_data = nullptr;
}

dso::HarmonicCoeffs &
dso::HarmonicCoeffs::operator=(dso::HarmonicCoeffs &&h) noexcept {
  this->m_degree = h.m_degree;
  this->m_data = h.m_data;
  h.m_degree = 0;
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

  return 0;
}
