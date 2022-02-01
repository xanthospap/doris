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
void print() noexcept {
    for (int i=0; i<=m_degree; i++) {
        for (int j=0; j<=m_degree; j++) {
            printf("%15.10e ", m_data[i][j]);
        }
        printf("\n");
    }
}
#endif

int degree() const noexcept {return m_degree;}

double *C_row(int degree) noexcept {
    return m_data[degree]; // C(degree,0)-> C(degree, degree)
}

double &C(int i, int j) noexcept {
    return C_row(i)[j];
}

//double &C2(int i, int j) noexcept {
//    int r = i;
//    int c = j;
//    printf("\tC(%02d, %02d) = data[%02d][%02d]\n", i, j, r, c);
//    return dummy;
//}

double *S_row(int degree) noexcept {
    int off = m_degree - degree;
    return m_data[off] + off + 1; // S(degree,1)-> S(degree, degree)
}

double &S(int i, int j) noexcept {
    return S_row(i)[j-1];
}

//double &S2(int i, int j) noexcept {
//    int off = m_degree - i;
//  int r = off;
//  int c = off + j;
//  printf("\tS(%02d, %02d) = data[%02d][%02d]\n", i, j, r, c);
//  return dummy;
//}

}; //HarmonicCoeffs

#endif
