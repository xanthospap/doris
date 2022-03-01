#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "atmosphere.hpp"
#include "doctest.h"
#include <cstdio>
#include <cstring>

using dso::air_density_models::msis86::msis86;

TEST_CASE("MSIS-86 Comparisson of results based on Fortran implementation.") {
  int day[] = {172, 81, 172, 172, 172, 172, 172, 172, 172, 172};
  double sec[] = {29000e0, 29000e0, 75000e0, 29000e0, 29000e0,
                  29000e0, 29000e0, 29000e0, 29000e0, 29000e0};
  double alt[] = {400e0, 400e0, 400e0, 100e0, 400e0,
                  400e0, 400e0, 400e0, 400e0, 400e0};
  double glat[] = {60e0, 60e0, 60e0, 60e0, 0e0, 60e0, 60e0, 60e0, 60e0, 60e0};
  double glong[] = {-70e0, -70e0, -70e0, -70e0, -70e0,
                  0e0,   -70e0, -70e0, -70e0, -70e0};
  double stl[] = {16e0, 16e0, 16e0, 16e0, 16e0, 16e0, 4e0, 16e0, 16e0, 16e0};
  double f107a[] = {150e0, 150e0, 150e0, 150e0, 150e0,
                  150e0, 150e0, 70e0,  150e0, 150e0};
  double f107[] = {150e0, 150e0, 150e0, 150e0, 150e0,
                 150e0, 150e0, 70e0,  180e0, 150e0};
  double ap[] = {4, 4, 4, 4, 4, 4, 4, 4, 4, 40e0, 0, 0, 0, 0, 0, 0};

  int switches[24];
  double outd[8], outt[2];
  std::memset(switches, 1, sizeof(int) * 24);

  const double results[][10] = {
      {1277.31, 1195.09, 1273.80, 1277.31, 1207.30, 1228.82, 1125.30, 1007.50,
       1331.63, 1384.07},
      {1265.63, 1188.85, 1262.36, 181.10, 1204.56, 1220.49, 1121.31, 1001.96,
       1316.23, 1367.98},
      {6519495.625, 33673157.500, 7184919.375, 1012165520.000, 18569566.250,
       8834627.500, 5582126.875, 3484452.812, 6649606.875, 5520741.875},
      {984292240.000, 1538254240.000, 1035035440.000, 3068923084800.000,
       1497053760.000, 1190952080.000, 610884320.000, 357530720.000,
       1095939600.000, 1045238560.000},
      {156377280.000, 143350910.000, 154863870.000, 82722611527680.000,
       166651870.000, 142276740.000, 95626520.000, 33023442.500, 188760840.000,
       304576960.000},
      {5490940.625, 5618148.125, 5428244.375, 19583061196800.000, 4545839.688,
       3966874.688, 3268814.375, 1054376.094, 6832877.500, 14832890.000},
      {32790.632, 19203.220, 29907.495, 926307942400.000, 17990.801, 22639.260,
       13863.903, 3404.718, 43365.361, 103536.016},
      {348313.242, 523080.781, 324106.191, 145393230.000, 615790.352,
       369626.602, 531904.219, 910571.250, 326974.746, 269174.219},
      {33079490.000, 39983952.500, 31845422.500, 4351691.875, 50382005.000,
       31282887.500, 8037972.500, 14231721.250, 45703915.000, 30068880.000},
      {0.00000000000003452, 0.00000000000004897, 0.00000000000003577,
       0.00000000502822339, 0.00000000000004905, 0.00000000000003924,
       0.00000000000002107, 0.00000000000001144, 0.00000000000003935,
       0.00000000000004345}};

  int i = 0;
  //SUBCASE("Test-Case 1") {
    msis86(day[i], sec[i], alt[i], glat[i], glong[i], stl[i], f107a[i], f107[i],
           switches, &ap[i], outd, outt);
    CHECK(std::abs(outt[0] - results[0][i]) < 1e-2);
    CHECK(std::abs(outt[1] - results[1][i]) < 1e-2);
    CHECK(std::abs(outd[0] - results[2][i]) < 1e-2);
    CHECK(std::abs(outd[1] - results[3][i]) < 1e-2);
    CHECK(std::abs(outd[2] - results[4][i]) < 1e-2);
    CHECK(std::abs(outd[3] - results[5][i]) < 1e-2);
    CHECK(std::abs(outd[4] - results[6][i]) < 1e-2);
    CHECK(std::abs(outd[5] - results[7][i]) < 1e-17);
    CHECK(std::abs(outd[6] - results[8][i]) < 1e-2);
    CHECK(std::abs(outd[7] - results[9][i]) < 1e-2);
    for (int k=0; k<8; k++) printf("d[%d] = %15.1f\n", k, outd[k]);
  //}
  
  i = 1;
  //SUBCASE("Test-Case 1") {
    msis86(day[i], sec[i], alt[i], glat[i], glong[i], stl[i], f107a[i], f107[i],
           switches, &ap[i], outd, outt);
    CHECK(std::abs(outt[0] - results[0][i]) < 1e-2);
    CHECK(std::abs(outt[1] - results[1][i]) < 1e-2);
    CHECK(std::abs(outd[0] - results[2][i]) < 1e-2);
    CHECK(std::abs(outd[1] - results[3][i]) < 1e-2);
    CHECK(std::abs(outd[2] - results[4][i]) < 1e-2);
    CHECK(std::abs(outd[3] - results[5][i]) < 1e-2);
    CHECK(std::abs(outd[4] - results[6][i]) < 1e-2);
    CHECK(std::abs(outd[5] - results[7][i]) < 1e-17);
    CHECK(std::abs(outd[6] - results[8][i]) < 1e-2);
    CHECK(std::abs(outd[7] - results[9][i]) < 1e-2);
  //}
  
  /*
  i = 2;
  SUBCASE("Test-Case 1") {
    msis86(day[i], sec[i], alt[i], glat[i], glong[i], stl[i], f107a[i], f107[i],
           switches, &ap[i], outd, outt);
    CHECK(std::abs(outt[0] - results[0][i]) < 1e-2);
    CHECK(std::abs(outt[1] - results[1][i]) < 1e-2);
    CHECK(std::abs(outd[0] - results[2][i]) < 1e-2);
    CHECK(std::abs(outd[1] - results[3][i]) < 1e-2);
    CHECK(std::abs(outd[2] - results[4][i]) < 1e-2);
    CHECK(std::abs(outd[3] - results[5][i]) < 1e-2);
    CHECK(std::abs(outd[4] - results[6][i]) < 1e-2);
    CHECK(std::abs(outd[5] - results[7][i]) < 1e-17);
    CHECK(std::abs(outd[6] - results[8][i]) < 1e-2);
    CHECK(std::abs(outd[7] - results[9][i]) < 1e-2);
  }
  
  i = 3;
  SUBCASE("Test-Case 1") {
    msis86(day[i], sec[i], alt[i], glat[i], glong[i], stl[i], f107a[i], f107[i],
           switches, &ap[i], outd, outt);
    CHECK(std::abs(outt[0] - results[0][i]) < 1e-2);
    CHECK(std::abs(outt[1] - results[1][i]) < 1e-2);
    CHECK(std::abs(outd[0] - results[2][i]) < 1e-2);
    CHECK(std::abs(outd[1] - results[3][i]) < 1e-2);
    CHECK(std::abs(outd[2] - results[4][i]) < 1e-2);
    CHECK(std::abs(outd[3] - results[5][i]) < 1e-2);
    CHECK(std::abs(outd[4] - results[6][i]) < 1e-2);
    CHECK(std::abs(outd[5] - results[7][i]) < 1e-17);
    CHECK(std::abs(outd[6] - results[8][i]) < 1e-2);
    CHECK(std::abs(outd[7] - results[9][i]) < 1e-2);
  }
  
  i = 4;
  SUBCASE("Test-Case 1") {
    msis86(day[i], sec[i], alt[i], glat[i], glong[i], stl[i], f107a[i], f107[i],
           switches, &ap[i], outd, outt);
    CHECK(std::abs(outt[0] - results[0][i]) < 1e-2);
    CHECK(std::abs(outt[1] - results[1][i]) < 1e-2);
    CHECK(std::abs(outd[0] - results[2][i]) < 1e-2);
    CHECK(std::abs(outd[1] - results[3][i]) < 1e-2);
    CHECK(std::abs(outd[2] - results[4][i]) < 1e-2);
    CHECK(std::abs(outd[3] - results[5][i]) < 1e-2);
    CHECK(std::abs(outd[4] - results[6][i]) < 1e-2);
    CHECK(std::abs(outd[5] - results[7][i]) < 1e-17);
    CHECK(std::abs(outd[6] - results[8][i]) < 1e-2);
    CHECK(std::abs(outd[7] - results[9][i]) < 1e-2);
  }
  
  i = 5;
  SUBCASE("Test-Case 1") {
    msis86(day[i], sec[i], alt[i], glat[i], glong[i], stl[i], f107a[i], f107[i],
           switches, &ap[i], outd, outt);
    CHECK(std::abs(outt[0] - results[0][i]) < 1e-2);
    CHECK(std::abs(outt[1] - results[1][i]) < 1e-2);
    CHECK(std::abs(outd[0] - results[2][i]) < 1e-2);
    CHECK(std::abs(outd[1] - results[3][i]) < 1e-2);
    CHECK(std::abs(outd[2] - results[4][i]) < 1e-2);
    CHECK(std::abs(outd[3] - results[5][i]) < 1e-2);
    CHECK(std::abs(outd[4] - results[6][i]) < 1e-2);
    CHECK(std::abs(outd[5] - results[7][i]) < 1e-17);
    CHECK(std::abs(outd[6] - results[8][i]) < 1e-2);
    CHECK(std::abs(outd[7] - results[9][i]) < 1e-2);
  }
  */
}