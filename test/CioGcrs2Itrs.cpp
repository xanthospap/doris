#include "iers2010/iau.hpp"
#include "orbit_integration.hpp"
#include <stdio.h>

int main() {
  const double tt1 = 2.45419550000000000e+06;
  const double tt2 = 5.00754444444444391e-01;
  const double ut11 = 2.45419550000000000e+06;
  const double ut12 = 4.99999165813830970e-01;
  const double dx06 = 8.48423941941688022e-10;
  const double dy06 = -1.09519410562644187e-09;
  const double xp = 1.69336692165300939e-07;
  const double yp = 2.34318354543240798e-06;

  double X, Y;
  iers2010::sofa::xy06(tt1, tt2, X, Y);
  const double s = iers2010::sofa::s06(tt1, tt2, X, Y);
  printf("// x = %+.9f\n", X);
  printf("// y = %+.9f\n", Y);
  printf("// s = %+.9f\n", s);
  X += dx06;
  Y += dy06;
  const Eigen::Matrix<double, 3, 3> rc2i = iers2010::sofa::c2ixys_e(X, Y, s);
  printf("// GCRS-to-CIRS");
  for (int i = 0; i < 3; i++) {
    printf("\n//");
    for (int j = 0; j < 3; j++) {
      printf(" %+.9f ", rc2i(i, j));
    }
  }
  const double era = iers2010::sofa::era00(ut11, ut12);
  printf("\n// era = %+.9f\n", era);
  const auto rc2ti = Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) * rc2i;
  printf("// GCRS-to-TIRS");
  for (int i = 0; i < 3; i++) {
    printf("\n//");
    for (int j = 0; j < 3; j++) {
      printf(" %+.9f ", rc2ti(i, j));
    }
  }
  const double sp = iers2010::sofa::sp00(tt1, tt2);
  const auto rpom = iers2010::sofa::pom00_e(xp, yp, sp);
  printf("\n// POM");
  for (int i = 0; i < 3; i++) {
    printf("\n//");
    for (int j = 0; j < 3; j++) {
      printf(" %+.9f ", rpom(i, j));
    }
  }

  const auto g2i = rpom * rc2ti;
  printf("\n// GCRS-to-ITRS");
  for (int i = 0; i < 3; i++) {
    printf("\n//");
    for (int j = 0; j < 3; j++) {
      printf(" %+.9f ", g2i(i, j));
    }
  }
  printf("\n");

  return 0;
}
