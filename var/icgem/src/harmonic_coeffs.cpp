#include "harmonic_coeffs.hpp"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cassert>

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

double fsum(int n, int m) noexcept {
  if (!m) return 1e0;
  int start = n-m+1;
  int stop = n+m;
  double sum = 1e0;
  for (int i=start; i<=stop; i++) sum *= (double)i;
  return sum;
}

int compute_fac(int n, int m, double *fac) noexcept {
  for (int i=0; i<=m; i++)
    fac[m] = std::sqrt(fsum(n,i) / (2*n+1) / (2-(i==0)));
  return 0;
}

int HarmonicCoeffs::denormalize(int order) noexcept {
  if (order < 0)
    order = m_degree;
  if (order > m_degree) {
    fprintf(stderr,
            "[ERROR] Invalid order provided (%d) for harmonics of max degree "
            "%d! Cannot denormalize (traceback: %s)\n",
            order, m_degree, __func__);
    return 1;
  }
  printf("Denormalizing with degree=%d and order=%d\n", m_degree, order);

  assert(m_degree>0);
  double *facs = new double[m_degree];
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
      s[m-1] /= facs[m-1];
  }
  delete[] facs;

  return 0;
}