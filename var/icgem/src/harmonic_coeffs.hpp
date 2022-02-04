#ifndef __HARMONIC_POTENTIAL_COEFFICIENTS_HPP__
#define __HARMONIC_POTENTIAL_COEFFICIENTS_HPP__

#ifdef DEBUG
#include <cstdio>
#endif

class HarmonicCoeffs {
public:
    int m_degree{0};
    double **m_data{nullptr};
    double dummy;

    int allocate() noexcept;
    int deallocate() noexcept;

public:
  HarmonicCoeffs() : m_degree(0), m_data(nullptr) {};
  HarmonicCoeffs(int n) : m_degree(n) { allocate(); }
  HarmonicCoeffs(const HarmonicCoeffs &h) = delete;
  HarmonicCoeffs &operator=(const HarmonicCoeffs &h) = delete;
  HarmonicCoeffs(HarmonicCoeffs &&h) noexcept;
  HarmonicCoeffs &operator=(HarmonicCoeffs &&h) noexcept;
  ~HarmonicCoeffs() noexcept { deallocate(); }

#ifdef DEBUG
void print(double scale=1e0) noexcept {
    for (int i=0; i<=m_degree; i++) {
        for (int j=0; j<=m_degree; j++) {
            printf("%15.10e ", scale * m_data[i][j]);
        }
        printf("\n");
    }
}
#endif

int degree() const noexcept {return m_degree;}

int denormalize(int order=-1) noexcept;

double *C_row(int degree) noexcept {
    return m_data[degree]; // C(degree,0)-> C(degree, degree)
}
const double *C_row(int degree) const noexcept {
    return m_data[degree]; // C(degree,0)-> C(degree, degree)
}

double &C(int i, int j) noexcept {
    return C_row(i)[j];
}

const double &C(int i, int j) const noexcept {
    return C_row(i)[j];
}

double *S_row(int degree) noexcept {
    int off = m_degree - degree;
    return m_data[off] + off + 1; // S(degree,1)-> S(degree, degree)
}

const double *S_row(int degree) const noexcept {
    int off = m_degree - degree;
    return m_data[off] + off + 1; // S(degree,1)-> S(degree, degree)
}

/// @warning never ask for a coefficient with order (aka j) = 0. All S_nm for
///          m=0 are equal to 0e0
double &S(int i, int j) noexcept {
    return S_row(i)[j-1];
}

/// @warning never ask for a coefficient with order (aka j) = 0. All S_nm for
///          m=0 are equal to 0e0
const double &S(int i, int j) const noexcept {
    return S_row(i)[j-1];
}

}; //HarmonicCoeffs

#endif
