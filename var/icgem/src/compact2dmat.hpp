#ifndef __COMAPCT_2D_SIMPLE_MATRIX_HPP__
#define __COMAPCT_2D_SIMPLE_MATRIX_HPP__

#include <cstring>
#include <algorithm>
#ifdef DEBUG
#include <cstdio>
#endif

enum class MatrixStorageType : char {RowWise, ColumnWise, Trapezoid};

template<MatrixStorageType S>
struct StorageImplementation {};

template <> struct StorageImplementation<MatrixStorageType::Trapezoid> {
  int rows, cols;
  constexpr StorageImplementation(int r, int c) noexcept : rows(r), cols(c){};

  constexpr std::size_t num_elements() const noexcept {
    if (rows == cols) {
      std::size_t n = rows;
      return (n*(n+1)) / 2;
    } else {
      std::size_t n = 0;
      for (int r=0; r<rows; r++)
        n += std::min(r+1, cols);
      return n;
    }
  }

  constexpr int pts_in_row(int row) const noexcept {
    return std::min(row+1, cols);
  }

  constexpr int slice(int row) const noexcept {
    int offset = 0;
    for (int i=0; i<row; i++)
      offset += pts_in_row(i);
    return offset;
  }
  constexpr int element_offset(int row, int column) const noexcept {
      return (column<=rows) ? (slice(row) + column) : (slice(column) + row);
  }
  #ifdef DEBUG
  void print(const double *data) const noexcept {
    for (int i=0; i<rows; i++) {
      for (int j=0; j<pts_in_row(i); j++) {
        printf("%15.10e ", data[element_offset(i,j)]);
      }
      printf("\n");
    }
    return;
  }
  #endif
};

template<>
struct StorageImplementation<MatrixStorageType::RowWise> {
    int rows, cols;
    constexpr StorageImplementation(int r, int c) noexcept: rows(r), cols(c) {};
    
    std::size_t num_elements() const noexcept {return rows * cols;}
    constexpr int  
    element_offset(int row, int column) const noexcept {
      return row * cols + column;
    }
    constexpr int slice(int row) const noexcept {
        return row * cols;
    }
  #ifdef DEBUG
  void print(const double *data) const noexcept {
    for (int i=0; i<rows; i++) {
      for (int j=0; j<cols; j++) {
        printf("%15.10e ", data[element_offset(i,j)]);
      }
      printf("\n");
    }
    return;
  }
  #endif
};

template<>
struct StorageImplementation<MatrixStorageType::ColumnWise> {
    int rows, cols;
    constexpr StorageImplementation(int r, int c) noexcept : rows(r), cols(c){};

    std::size_t num_elements() const noexcept {return rows * cols;}
    constexpr int element_offset(int row, int column) const noexcept {
      return column * rows + row;
    }
    constexpr int slice(int col) const noexcept {
        return col * rows;
    }
  #ifdef DEBUG
  void print(const double *data) const noexcept {
    for (int i=0; i<rows; i++) {
      for (int j=0; j<cols; j++) {
        printf("%15.10e ", data[element_offset(i,j)]);
      }
      printf("\n");
    }
    return;
  }
  #endif
};

template<MatrixStorageType S>
class Mat2D {
private:
  StorageImplementation<S> m_storage;
  double *m_data{nullptr};

public :
  int rows() const noexcept {return m_storage.rows;}
  int cols() const noexcept {return m_storage.cols;}
  long num_elements() const noexcept {return m_storage.num_elements();}
  double &operator()(int i, int j) noexcept { return m_data[m_storage.element_offset(i,j)]; }
  const double &operator()(int i, int j) const  noexcept { return m_data[m_storage.element_offset(i,j)]; }
  const double *slice(int i) const noexcept { return m_data[m_storage(slice(i))]; }
  double *slice(int i) noexcept { return m_data[m_storage(slice(i))]; }
  void fill_with(double val) noexcept {
    std::fill(m_data, m_data+m_storage.num_elements(), val);
  }
  const double *data() const noexcept {return m_data;}
  #ifdef DEBUG
  double *data() noexcept {return m_data;}
  #endif

  Mat2D(int rows, int cols) noexcept
      : m_storage(rows, cols), m_data(new double[m_storage.num_elements()]){};
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
  #ifdef DEBUG
  void print() const noexcept { return m_storage.print(m_data); }
  #endif

}; // Mat2D
  

#endif