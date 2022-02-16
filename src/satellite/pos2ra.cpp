#include "satellite.hpp"
#include <cmath>
#ifdef ASSERT_ERROR
#include <cassert>
#endif

int dso::pos2ad(const Vector3 &pos, double &a, double &d) noexcept {
  // calculate the magnitude of r
  const double r = pos.norm();
  
  // Calculate the direction cosines of r
  const double l = pos.x() / r;
  const double m = pos.y() / r;
  const double n = pos.z() / r;
  
  // Calculate the declination
#ifdef ASSERT_ERROR
  assert(n>=-1e0 && n<=1e0);
#endif
  d = std::asin(n); // range [-pi/2 to pi/2]
  
  // Calculate the right ascension
  const double af = std::acos(l / std::cos(d));
#ifdef BRANCHLESS
  const double ret[] = {iers2010::D2PI - af, af};
  a = ret[m>0];
#else
  a = (m<=0) ? (iers2010::D2PI - af) : af;
#endif

  return 0;
}
