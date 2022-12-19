#include "associated_legendre.hpp"
#include <cmath>

// See 
// To be used with the standard forward column method(s)
void dso::AssociatedLegendreFunctions::compute_factors() noexcept {
  const int maxDegree = m_degree;
  
  // diagonal terms (Eq. 13)
  F1->operator()(1,1) = std::sqrt(3e0);
  for (int n=2; n<=maxDegree; n++) {
    (*F1)(n,n) = std::sqrt( (2*n+1) / (2e0*n) );
  }

  // degrees > 0 and order m=0
  for (int n = 1; n <= maxDegree; n++) {
    const long ln = n;
    const long Anm_u = (2L * ln - 1L) * (2L * ln + 1L);
    const long Anm_d = ln * ln;
    const long Bnm_u = (2L * ln + 1L) * (ln - 1L) * (ln - 1L);
    const long Bnm_d = ln * ln * (2l * ln - 3);
    F1->operator()(n, 0) = std::sqrt((double)Anm_u / (double)Anm_d);
    F2->operator()(n, 0) = std::sqrt((double)Bnm_u / (double)Bnm_d);
  }

  // all other terms, per order (m), Eq. 12 (a_nm and b_nm)
  for (int m=1; m<=maxDegree; m++) {
    for (int n=m+1; n<=maxDegree; n++) {
      const long ln = n;
      const long lm = m;
      const long Anm_u = (2L*ln-1L)*(2L*ln+1L);
      const long Anm_d = (ln-lm)*(ln+lm);
      const long Bnm_u = (2L*ln+1L)*(ln+lm-1L)*(ln-lm-1L);
      const long Bnm_d = (ln-lm)*(ln+lm)*(2l*ln-3);
      F1->operator()(n,m) = std::sqrt((double)Anm_u/(double)Anm_d);
      F2->operator()(n,m) = std::sqrt((double)Bnm_u/(double)Bnm_d);

      //const long ipar = (n + m) * (n - m);
      //const double f = (2 * n + 1) / static_cast<double>(ipar);
      //F1->operator()(n,m) = std::sqrt(f*(2*n-1));
      //F2->operator()(n,m) = std::sqrt(f*((n-m-1)*(n+m-1))/(2*n-3));
    }
  }

  return;
}

/*void dso::AssociatedLegendreFunctions::compute(double angle) noexcept {
  const int maxDegree = m_degree;
  const double ca = std::cos(angle);
  const double sa = std::sin(angle);
  
  P(0,0) = 1e0;
  P(1,0) = 
}*/

// Standard Forward Column Method
void dso::AssociatedLegendreFunctions::compute(double angle) noexcept {
  const int maxDegree = m_degree;
  const double u = std::sin(angle);
  const double t = std::cos(angle);
  
  P(0,0) = 1e0;
  
  { // fill elements for order = 0
    const int m = 0;
    P(1,0) = (*F1)(1,0) * t * P(0,0);
    for (int n=2; n<=maxDegree; n++)
      P(n,m) = (*F1)(n,m) * t * P(n-1, m) - (*F2)(n,m) * P(n-2, m);
  }

  // column wise visiting, for orders: 1<=m<maxDegree
  double lastDiagonal = P(0,0);
  for (int m=1; m<maxDegree; m++) {
    // diagonal elements
    P(m,m) = (*F1)(m,m) * u * /*P(m-1,m-1)*/lastDiagonal;
    lastDiagonal = P(m,m);
    // sub-diagonal points
    P(m+1,m) = (*F1)(m+1,m) * t * P(m,m);
    // all other points of same order = m
    for (int n=m+2; n<=maxDegree; n++)
      P(n,m) = (*F1)(n,m) * t * P(n-1, m) - (*F2)(n,m) * P(n-2, m);
  }

  // we are missing the last (diagonal) element!
  {
    const int m = maxDegree;
    P(m, m) = (*F1)(m,m) * u * /*P(m-1,m-1)*/lastDiagonal;
  }

  // all done!
  return;
}
