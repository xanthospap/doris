#include <cassert>
#include <cmath>

/// Compute renormalized associated Legendre polynomials. m and l are integers 
/// satisfying 0 <= m <= l, while x is in [-1, 1]
double plegendre(int l, int m, double x) noexcept {
  assert(m >=0 && m <= l);
  assert(std::abs(x) <= 1e0);

  double omx2, fact;
  double pmm=1e0;
  if (m > 0) {
    omx2=(1.0-x)*(1.0+x);
    fact=1e0;
    for (int i=1;i<=m;i++) {
      pmm *= omx2*fact/(fact+1e0);
      fact += 2e0;
    }
  }

  pmm=std::sqrt((2e0*m+1e0)*pmm/(4e0*ngpt::DPI));
  if (m & 1)
    pmm=-pmm;

  if (l == m)
    return pmm;
  else {
    double pmmp1=x*std::sqrt(2e0*m+3e0)*pmm;
    if (l == (m+1))
      return pmmp1;
    else {
        double oldfact=std::sqrt(2e0*m+3e0);
        double pll;
        for (int ll=m+2;ll<=l;ll++) {
          fact=std::sqrt((4e0*ll*ll-1e0)/(ll*ll-m*m));
          pll=(x*pmmp1-pmm/oldfact)*fact;
          oldfact=fact;
          pmm=pmmp1;
          pmmp1=pll;
        }
        return pll;
    }
  }
  
}
