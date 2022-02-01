#include "compact2dmat.hpp"
#include <cstring>

void Mat2D::copy_data(const double *src, long N) noexcept {
    std::memcpy(this->m_data, src, sizeof(double)*N);
}

Mat2D &Mat2D::operator=(const Mat2D &mat) noexcept {
  if (this != &mat) {
    long N = mat.m_rows * mat.m_cols;
    if (m_data)
      delete[] m_data;
    m_data = new double[N];
    if (m_data) {
      copy_data(mat.m_data, N);
      m_rows = mat.m_rows;
      m_cols = mat.m_cols;
    } else {
      m_rows = 0;
      m_cols = 0;
    }
  }
  return *this;
}

Mat2D &Mat2D::operator=(Mat2D &&mat) noexcept {
    if (m_data)
      delete[] m_data;
    
    m_data = mat.m_data;
    m_rows = mat.m_rows;
    m_cols = mat.m_cols;

    mat.m_data = nullptr;
    mat.m_rows = 0;
    mat.m_cols = 0;

  return *this;
}