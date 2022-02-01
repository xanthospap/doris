#ifndef __COMAPCT_2D_SIMPLE_MATRIX_HPP__
#define __COMAPCT_2D_SIMPLE_MATRIX_HPP__

enum class MatrixStorageType : char {RowWise, ColumnWise, Triangular};

template<MatrixStorageType S>
struct StorageImplementation {};

template<>
struct StorageImplementation<MatrixStorageType::RowWise> {
    int rows, cols;
    constexpr long element_offset(int row, int column) const noexcept {
        return row * cols + column;
    }
    constexpr int row_offset(int row) const noexcept {
        return row * cols;
    }
};

template<>
struct StorageImplementation<MatrixStorageType::ColumnWise> {
    int rows, cols;
    constexpr long element_offset(int row, int column) const noexcept {
        return row * cols + column;
    }
    constexpr int col_offset(int row) const noexcept {
        return row * cols;
    }
};

template<MatrixStorageType S>
class Mat2D {
private:
  int m_rows{0};
  int m_cols{0};
  double *m_data{nullptr};

  void copy_data(const double *src, long size) noexcept;

  Mat2D(int rows, int cols) noexcept
      : m_rows(rows), m_cols(cols), m_data(new double[rows * cols]){};
  ~Mat2D() noexcept {
    if (m_data)
      delete[] m_data;
  };
  Mat2D(const Mat2D &m) noexcept
      : m_rows(m.m_rows), m_cols(m.m_cols),
        m_data(new double[m.m_rows * m.m_cols]) {
    copy_data(m.m_data, m.m_rows * m.m_cols);
  };
  Mat2D(Mat2D &&m) noexcept
      : m_rows(m.m_rows), m_cols(m.m_cols), m_data(m.m_data) {
    m.m_rows = 0;
    m.m_cols = 0;
    m.m_data = nullptr;
  };
  Mat2D &operator=(const Mat2D &mat) noexcept;
  Mat2D &operator=(Mat2D &&mat) noexcept;

public:
  double *row(int i) noexcept { return m_data + i * m_cols; }
  const double *row(int i) const noexcept { return m_data + i * m_cols; }
  double &operator()(int i, int j) noexcept {
    return *(m_data + i * m_cols + j);
  }
}; // Mat2D

#endif