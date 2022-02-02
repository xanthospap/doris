#ifndef __COMAPCT_2D_SIMPLE_MATRIX_HPP__
#define __COMAPCT_2D_SIMPLE_MATRIX_HPP__

#include <cstring>

enum class MatrixStorageType : char {RowWise, ColumnWise, Triangular};

template<MatrixStorageType S>
struct StorageImplementation {};

template<>
struct StorageImplementation<MatrixStorageType::RowWise> {
    int rows, cols;
    constexpr StorageImplementation(int r, int c) noexcept: rows(r), cols(c) {};
    
    constexpr int  
    element_offset(int row, int column) const noexcept {
      return row * cols + column;
    }
    constexpr int slice(int row) const noexcept {
        return row * cols;
    }
};

template<>
struct StorageImplementation<MatrixStorageType::ColumnWise> {
    int rows, cols;
    constexpr StorageImplementation(int r, int c) noexcept : rows(r), cols(c){};

    constexpr int element_offset(int row, int column) const noexcept {
      return column * rows + row;
    }
    constexpr int slice(int col) const noexcept {
        return col * rows;
    }
};

template<MatrixStorageType S>
class Mat2D {
private:
  StorageImplementation<S> m_storage;
  double *m_data{nullptr};

public :
  int rows() const noexcept {return m_storage.rows;}
  int cols() const noexcept {return m_storage.cols;}
  long num_elements() const noexcept {return m_storage.rows * m_storage.cols;}
  double &operator()(int i, int j) noexcept { return m_data[m_storage.element_offset(i,j)]; }
  const double &operator()(int i, int j) const  noexcept { return m_data[m_storage.element_offset(i,j)]; }
  const double *slice(int i) const noexcept { return m_data[m_storage(slice(i))]; }
  double *slice(int i) noexcept { return m_data[m_storage(slice(i))]; }

  Mat2D(int rows, int cols) noexcept
      : m_storage(rows, cols), m_data(new double[rows * cols]){};
  ~Mat2D() noexcept {
    if (m_data)
      delete[] m_data;
  }
  Mat2D(const Mat2D &mat) noexcept
      : m_storage(mat.m_storage), m_data(new double[mat.num_elements()]) {
    std::memset(m_data, mat.m_data, sizeof(double) * mat.num_elements());
  }
  Mat2D(Mat2D &&mat) noexcept : m_storage(mat.m_storage), m_data(mat.m_data) {
    mat.m_data = nullptr;
  }
  Mat2D &operator=(const Mat2D &mat) noexcept {
    if (this!=&mat) {
      if (m_data && (m_storage.num_elements() != mat.num_elements())) {
        delete[] m_data;
        m_data = new double[mat.num_elements()];
      }
      std::memset(m_data, mat.m_data, sizeof(double) * mat.num_elements());
      m_storage.rows = mat.rows();
      m_storage.cols = mat.cols();
    }
    return *this;
  }
  Mat2D &operator=(Mat2D &&mat) noexcept { 
    m_storage.rows = mat.rows();
    m_storage.cols = mat.cols();
    m_data = mat.m_data;
    mat.m_data = nullptr;
  }

}; // Mat2D

#endif