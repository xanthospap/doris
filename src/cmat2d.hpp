#ifndef __COMPACT_2D_SIMPLE_MATRIX_HPP__
#define __COMPACT_2D_SIMPLE_MATRIX_HPP__

#include <algorithm>
#include <cstring>
#ifdef DEBUG
#include <cassert>
#include <cstdio>
#endif

/// @file Naive implementation of 2-Dimensional matrices; note that this
/// implementation targets 2-dimensional matrices that are only meant to store
/// data NOT perform arithmetic operations.

namespace dso {

///< Enum class to describe storage type of a 2-d matrix
enum class MatrixStorageType : char {
  RowWise,             ///< Row-Wise storage
  ColumnWise,          ///< Column-Wise storage
  Trapezoid,           ///< A trapezoid matrix with row-wise storage
  LwTriangularRowWise, ///< Lower triangular, stored Row-Wise
  LwTriangularColWise  ///< Lower triangular, stored Col-Wise
};                     // MatrixStorageType

/// @brief implementation details depending on storage type, aka
///        MatrixStorageType
template <MatrixStorageType S> struct StorageImplementation {};

/// @brief Implementation details for a 2-d lower triangular matrix, holding
///        data in a Row-Wise fashion.
/// Here we are not interested on the actual data of the matrix, but only the
/// indexing implementation of the matrix.
template <>
struct StorageImplementation<MatrixStorageType::LwTriangularRowWise> {
  int rows;
  constexpr StorageImplementation(int r, [[maybe_unused]] int _) noexcept
      : rows(r){};

  /// @brief Compute number of elements stored
  constexpr std::size_t num_elements() const noexcept {
    return rows * (rows + 1) / 2;
  }

  constexpr std::size_t _size(int d) const noexcept { return d * (d + 1) / 2; }

  constexpr int nrows() const noexcept { return rows; }
  constexpr int ncols() const noexcept { return rows; }

  /// Return the offset from the begining of the data array, given a row
  /// number. First row is row 0 (NOT row 1).
  /// That means that if the data is stored in an array e.g.
  ///   double *data = new double[num_pts];
  ///   double *row_3 = data[0] + slice(2);
  /// will point to the first (0) element of the third row.
  constexpr int slice(int row) const noexcept {
    return _size(row - 1) * (row != 0);
  }

  /// @brief Index of element (row, column) in the data array.
  /// E.g. data[element_offset(1,2)] will return the element in the second
  /// row, and third column.
  constexpr int element_offset(int row, int column) const noexcept {
    return slice(row) + column;
  }
}; // StorageImplementation<MatrixStorageType::LwTriangularRowWise>

/// @brief Implementation details for a 2-d lower triangular matrix, holding
///        data in a Col-Wise fashion.
/// Here we are not interested on the actual data of the matrix, but only the
/// indexing implementation of the matrix.
template <>
struct StorageImplementation<MatrixStorageType::LwTriangularColWise> {
  int rows;
  constexpr StorageImplementation(int r, [[maybe_unused]] int _) noexcept
      : rows(r){};

  /// @brief Compute number of elements stored
  constexpr std::size_t num_elements() const noexcept {
    return rows * (rows + 1) / 2;
  }

  constexpr int nrows() const noexcept { return rows; }
  constexpr int ncols() const noexcept { return rows; }

  /// @brief Index/offset of given column.
  /// Return the offset from the begining of the data array, given a column
  /// number.
  /// First column is column 0 (NOT row 1).
  /// That means that if the data is stored in an array e.g.
  ///   double *data = new double[num_pts];
  ///   double *col_3 = data[0] + slice(2);
  /// will point to the first (0) element of the third column.
  constexpr int slice(int col) const noexcept {
    return col * rows - col * (col - 1) / 2;
  }

  /// @brief Index of element (row, column) in the data array.
  /// E.g. data[element_offset(1,2)] will return the element in the second
  /// row, and third column.
  constexpr int element_offset(int row, int col) const noexcept {
    return slice(col) + (row - col);
  }
}; // StorageImplementation<MatrixStorageType::LwTriangularColWise>

/// @brief Implementation details for a 2-d trapezoid matrix, holding data in
///        a Row-Wise fashion.
/// In case rows == columns, this is actually a triangular matrix.
/// Here we are not interested on the actual data of the matrix, but only the
/// indexing implementation of the matrix.
template <> struct StorageImplementation<MatrixStorageType::Trapezoid> {
  int rows, cols;
  constexpr StorageImplementation(int r, int c) noexcept : rows(r), cols(c){};

  constexpr int nrows() const noexcept { return rows; }
  constexpr int ncols() const noexcept { return cols; }

  /// @brief Compute number of elements stored
  constexpr std::size_t num_elements() const noexcept {
    if (rows == cols) {
      std::size_t n = rows;
      return (n * (n + 1)) / 2;
    } else {
      std::size_t n = 0;
      for (int r = 0; r < rows; r++)
        n += std::min(r + 1, cols);
      return n;
    }
  }

  /// @brief Number of elements (data points) stored for a given row
  constexpr int pts_in_row(int row) const noexcept {
    return std::min(row + 1, cols);
  }

  /// Return the offset from the begining of the data array, given a row
  /// number. First row is row 0 (NOT row 1).
  /// That means that if the data is stored in an array e.g.
  ///   double *data = new double[num_pts];
  ///   double *row_3 = data[0] + slice(2);
  /// will point to the first (0) element of the third row.
  constexpr int slice(int row) const noexcept {
    int offset = 0;
    for (int i = 0; i < row; i++)
      offset += pts_in_row(i);
    return offset;
  }

  /// @brief Index of element (row, column) in the data array.
  /// E.g. data[element_offset(1,2)] will return the element in the second row,
  /// and third column.
  constexpr int element_offset(int row, int column) const noexcept {
    return (column <= rows) ? (slice(row) + column) : (slice(column) + row);
  }

}; // StorageImplementation<MatrixStorageType::Trapezoid>

/// @brief Implementation details for a 2-d dense matrix, holding data in
///        a Row-Wise fashion.
/// Here we are not interested on the actual data of the matrix, but only the
/// indexing implementation of the matrix.
template <> struct StorageImplementation<MatrixStorageType::RowWise> {
  int rows, cols;
  constexpr StorageImplementation(int r, int c) noexcept : rows(r), cols(c){};

  constexpr int nrows() const noexcept { return rows; }
  constexpr int ncols() const noexcept { return rows; }

  /// @brief Number of elements in matrix
  constexpr std::size_t num_elements() const noexcept { return rows * cols; }

  /// @brief Index of element (row, column) in the data array.
  /// E.g. data[element_offset(1,2)] will return the element in the second
  /// row, and third column.
  constexpr int element_offset(int row, int column) const noexcept {
    return row * cols + column;
  }

  /// @brief Index/offset of given row.
  /// Return the offset from the begining of the data array, given a row
  /// number.
  /// First row is row 0 (NOT row 1).
  /// That means that if the data is stored in an array e.g.
  ///   double *data = new double[num_pts];
  ///   double *row_3 = data[0] + slice(2);
  /// will point to the first (0) element of the third row.
  constexpr int slice(int row) const noexcept { return row * cols; }
}; // StorageImplementation<MatrixStorageType::RowWise>

/// @brief Implementation details for a 2-d dense matrix, holding data in
///        a Column-Wise fashion.
/// Here we are not interested on the actual data of the matrix, but only the
/// indexing implementation of the matrix.
template <> struct StorageImplementation<MatrixStorageType::ColumnWise> {
  int rows, cols;
  constexpr StorageImplementation(int r, int c) noexcept : rows(r), cols(c){};

  constexpr int nrows() const noexcept { return rows; }
  constexpr int ncols() const noexcept { return rows; }

  /// @brief Number of elements in matrix
  constexpr std::size_t num_elements() const noexcept { return rows * cols; }

  /// @brief Index of element (row, column) in the data array.
  /// E.g. data[element_offset(1,2)] will return the element in the second
  /// row, and third column.
  constexpr int element_offset(int row, int column) const noexcept {
    return column * rows + row;
  }

  /// @brief Index/offset of given column.
  /// Return the offset from the begining of the data array, given a column
  /// number.
  /// First column is column 0 (NOT row 1).
  /// That means that if the data is stored in an array e.g.
  ///   double *data = new double[num_pts];
  ///   double *col_3 = data[0] + slice(2);
  /// will point to the first (0) element of the third column.
  constexpr int slice(int col) const noexcept { return col * rows; }

}; // StorageImplementation<MatrixStorageType::ColumnWise>

/// @brief A naive implementation of a 2-d dense matrix.
/// The main objective of this class is to store data (e.g. coefficients) and
/// not perform arithmetic operations.
template <MatrixStorageType S> class Mat2D {
private:
  StorageImplementation<S> m_storage; ///< storage type; dictates indexing
  double *m_data{nullptr};            ///< the actual data

public:
  constexpr int rows() const noexcept { return m_storage.nrows(); }

  constexpr int cols() const noexcept { return m_storage.ncols(); }

  constexpr long num_elements() const noexcept {
    return m_storage.num_elements();
  }

  /// @brief Element indexing (rows and columns start from 0 --not 1--)
  double &operator()(int i, int j) noexcept {
#ifdef DEBUG
    assert(i < rows() && j < cols());
    assert(m_storage.element_offset(i, j) >= 0 &&
           m_storage.element_offset(i, j) < (int)m_storage.num_elements());
#endif
    return m_data[m_storage.element_offset(i, j)];
  }

  /// @brief Element indexing (rows and columns start from 0 --not 1--)
  const double &operator()(int i, int j) const noexcept {
#ifdef DEBUG
    assert(i < rows() && j < cols());
    assert(m_storage.element_offset(i, j) >= 0 &&
           m_storage.element_offset(i, j) < (int)m_storage.num_elements());
#endif
    return m_data[m_storage.element_offset(i, j)];
  }

  /// @brief Row/Column indexing (rows and columns start from 0 --not 1--)
  /// If the data is tored in a Row-Wise manner, this function will return the
  /// offset of the ith row; if data is stored in a column-wise manner, it will
  /// return the index of the first elelement of the ith column
  const double *slice(int i) const noexcept {
    return m_data + m_storage.slice(i);
  }

  /// @brief Row/Column indexing (rows and columns start from 0 --not 1--)
  /// If the data is tored in a Row-Wise manner, this function will return the
  /// offset of the ith row; if data is stored in a column-wise manner, it will
  /// return the index of the first elelement of the ith column
  double *slice(int i) noexcept { return m_data + m_storage.slice(i); }

  void fill_with(double val) noexcept {
    std::fill(m_data, m_data + m_storage.num_elements(), val);
  }

  void multiply(double value) noexcept {
    std::transform(m_data, m_data + m_storage.num_elements(), m_data,
                   [=](double d) { return d * value; });
  }

  const double *data() const noexcept { return m_data; }

#ifdef DEBUG
  double *data() noexcept { return m_data; }
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
    if (this != &mat) {
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

} // namespace dso

#endif
