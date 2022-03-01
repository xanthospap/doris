#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "atmosphere.hpp"
#include "doctest.h"
#include <cstdio>
#include <cstring>

using dso::air_density_models::msis86::msis86;

TEST_CASE("MSIS-86 Comparisson of results based on Matlab implementation.") {
  int day = 172;       // day of the year
  double sec = 29000.; // utc in seconds
  double alt = 400.;   // altitude in km
  double glat = 60.;   // latitude geodetic in degrees
  double glong = -70.; // longitude in degrees
  double stl = 16.;    // local solar time in hours
  double f107a = 150.; // mean solar flux over 90 days
  double f107 = 150.;  // solar flux from previous day
  double ap[] = {4, 4, 4, 4, 4, 4, 4}; // see note below
  
  int switches[24];
  double outd[8], outt[2];
  std::memset(switches, 1, sizeof(int) * 24);
  
  // note on ap array:
  // dayly Ap
  // 3 hr ap index for current time
  // 3 hr ap index for 3 hrs before current time
  // 3 hr ap index for 6 hrs before current time
  // 3 hr ap index for 9 hrs before current time
  // average of eight 3 hr ap indicies from 12 to 33 hrs prior to current time
  // average of eight 3 hr ap indicies from 36 to 59 hrs prior to current time

  SUBCASE("Example 1, testing for precesion 1e-3.") {
    msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
           outt);
    REQUIRE(std::abs(outd[0] - 667590.662) < 1e-3);
    REQUIRE(std::abs(outd[1] - 108792495.500) < 1e-3);
    REQUIRE(std::abs(outd[2] - 20702571.921) < 1e-3);
    REQUIRE(std::abs(outd[3] - 672686.058) < 1e-3);
    REQUIRE(std::abs(outd[4] - 4229.982) < 1e-3);
    REQUIRE(std::abs(outd[5] - 3.976e-15) < 1e-17);
    REQUIRE(std::abs(outd[6] - 34983.050) < 1e-3);
    REQUIRE(std::abs(outd[7] - 3609134.381) < 1e-3);
    REQUIRE(std::abs(outt[0] - 1277.3131) < 1e-3);
    REQUIRE(std::abs(outt[1] - 1270.0803) < 1e-3);
  }

  day = 81;
  SUBCASE("Example 2, testing for precesion 1e-3.") {
    msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
           outt);
    REQUIRE(std::abs(outd[0] - 3364853.653) < 1e-3);
    REQUIRE(std::abs(outd[1] - 153355366.963) < 1e-3);
    REQUIRE(std::abs(outd[2] - 14212339.605) < 1e-3);
    REQUIRE(std::abs(outd[3] - 558349.879) < 1e-3);
    REQUIRE(std::abs(outd[4] - 1905.489) < 1e-3);
    REQUIRE(std::abs(outd[5] - 4.88e-15) < 1e-17);
    REQUIRE(std::abs(outd[6] - 52300.580) < 1e-3);
    REQUIRE(std::abs(outd[7] - 3987746.529) < 1e-3);
    REQUIRE(std::abs(outt[0] - 1195.0909) < 1e-3);
    REQUIRE(std::abs(outt[1] - 1188.7496) < 1e-3);
  }

  day = 172;
  sec = 75000.;
  SUBCASE("Example 3, testing for precesion 1e-3.") {
    msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
           outt);
    REQUIRE(std::abs(outd[0] - 735697.269) < 1e-3);
    REQUIRE(std::abs(outd[1] - 114372564.992) < 1e-3);
    REQUIRE(std::abs(outd[2] - 20492466.708) < 1e-3);
    REQUIRE(std::abs(outd[3] - 664638.686) < 1e-3);
    REQUIRE(std::abs(outd[4] - 3855.342) < 1e-3);
    REQUIRE(std::abs(outd[5] - 4.11e-15) < 1e-17);
    REQUIRE(std::abs(outd[6] - 32552.226) < 1e-3);
    REQUIRE(std::abs(outd[7] - 3473758.383) < 1e-3);
    REQUIRE(std::abs(outt[0] - 1273.7997) < 1e-3);
    REQUIRE(std::abs(outt[1] - 1266.7298) < 1e-3);
  }

  sec = 29000.;
  alt = 100.;
  ap[6] = 40;
  SUBCASE("Example 4, testing for precesion 1e-3.") {
    msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
           outt);
    REQUIRE(std::abs(outd[0] - 103548323.467) < 1e-3);
    REQUIRE(std::abs(outd[1] - 304600073728.932) < 1e-3);
    REQUIRE(std::abs(outd[2] - 8465894630345.253) < 1e-3);
    REQUIRE(std::abs(outd[3] - 1973033277954.449) < 1e-3);
    REQUIRE(std::abs(outd[4] - 91901088091.951) < 1e-3);
    REQUIRE(std::abs(outd[5] - 5.12e-10) < 1e-17);
    REQUIRE(std::abs(outd[6] - 15250346.097) < 1e-3);
    REQUIRE(std::abs(outd[7] - 429826.936) < 1e-3);
    REQUIRE(std::abs(outt[0] - 1277.3131) < 1e-3);
    REQUIRE(std::abs(outt[1] - 180.2007) < 1e-3);
  }

  alt = 400.;
  glat = 0.;
  ap[5] = 40;
  ap[6] = 0;
  SUBCASE("Example 5, testing for precesion 1e-3.") {
    msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
           outt);
    REQUIRE(std::abs(outd[0] - 1856958.465) < 1e-3);
    REQUIRE(std::abs(outd[1] - 149705464.596) < 1e-3);
    REQUIRE(std::abs(outd[2] - 16665159.655) < 1e-3);
    REQUIRE(std::abs(outd[3] - 454583.228) < 1e-3);
    REQUIRE(std::abs(outd[4] - 1799.078) < 1e-3);
    REQUIRE(std::abs(outd[5] - 4.90e-15) < 1e-17);
    REQUIRE(std::abs(outd[6] - 61579.062) < 1e-3);
    REQUIRE(std::abs(outd[7] - 5038188.536) < 1e-3);
    REQUIRE(std::abs(outt[0] - 1207.2972) < 1e-3);
    REQUIRE(std::abs(outt[1] - 1204.5592) < 1e-3);
  }

  glat = 60.;
  glong = 0.;
  ap[4] = 40;
  ap[5] = 0;
  SUBCASE("Example 6, testing for precesion 1e-3.") {
    msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
           outt);
    REQUIRE(std::abs(outd[0] - 892286.980) < 1e-3);
    REQUIRE(std::abs(outd[1] - 126946511.216) < 1e-3);
    REQUIRE(std::abs(outd[2] - 17854729.479) < 1e-3);
    REQUIRE(std::abs(outd[3] - 458032.071) < 1e-3);
    REQUIRE(std::abs(outd[4] - 2721.189) < 1e-3);
    REQUIRE(std::abs(outd[5] - 4.31e-15) < 1e-17);
    REQUIRE(std::abs(outd[6] - 36773.619) < 1e-3);
    REQUIRE(std::abs(outd[7] - 3301014.134) < 1e-3);
    REQUIRE(std::abs(outt[0] - 1228.8197) < 1e-3);
    REQUIRE(std::abs(outt[1] - 1223.6185) < 1e-3);
  }

  glong = -70.;
  stl = 4.;
  ap[3] = 40;
  ap[4] = 0;
  SUBCASE("Example 7, testing for precesion 1e-3.") {
    msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
           outt);
    REQUIRE(std::abs(outd[0] - 562899.933) < 1e-3);
    REQUIRE(std::abs(outd[1] - 64622492.023) < 1e-3);
    REQUIRE(std::abs(outd[2] - 11833349.336) < 1e-3);
    REQUIRE(std::abs(outd[3] - 371380.789) < 1e-3);
    REQUIRE(std::abs(outd[4] - 1632.704) < 1e-3);
    REQUIRE(std::abs(outd[5] - 2.31e-15) < 1e-17);
    REQUIRE(std::abs(outd[6] - 52919.728) < 1e-3);
    REQUIRE(std::abs(outd[7] - 842655.227) < 1e-3);
    REQUIRE(std::abs(outt[0] - 1125.3016) < 1e-3);
    REQUIRE(std::abs(outt[1] - 1122.9512) < 1e-3);
  }

  stl = 16.;
  f107a = 70.;
  ap[2] = 40;
  ap[3] = 0;
  SUBCASE("Example 8, testing for precesion 1e-3.") {
    msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
           outt);
    REQUIRE(std::abs(outd[0] - 353888.980) < 1e-3);
    REQUIRE(std::abs(outd[1] - 39005044.366) < 1e-3);
    REQUIRE(std::abs(outd[2] - 4318005.925) < 1e-3);
    REQUIRE(std::abs(outd[3] - 127607.063) < 1e-3);
    REQUIRE(std::abs(outd[4] - 434.099) < 1e-3);
    REQUIRE(std::abs(outd[5] - 1.28e-15) < 1e-17);
    REQUIRE(std::abs(outd[6] - 90690.229) < 1e-3);
    REQUIRE(std::abs(outd[7] - 1532442.466) < 1e-3);
    REQUIRE(std::abs(outt[0] - 1007.5024) < 1e-3);
    REQUIRE(std::abs(outt[1] - 1004.5139) < 1e-3);
  }

  f107a = 150.;
  f107 = 180.;
  ap[1] = 40;
  ap[2] = 0;
  SUBCASE("Example 9, testing for precesion 1e-3.") {
    msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
           outt);
    REQUIRE(std::abs(outd[0] - 672454.105) < 1e-3);
    REQUIRE(std::abs(outd[1] - 117656963.133) < 1e-3);
    REQUIRE(std::abs(outd[2] - 24011903.418) < 1e-3);
    REQUIRE(std::abs(outd[3] - 801455.261) < 1e-3);
    REQUIRE(std::abs(outd[4] - 5317.831) < 1e-3);
    REQUIRE(std::abs(outd[5] - 4.40e-15) < 1e-17);
    REQUIRE(std::abs(outd[6] - 32519.041) < 1e-3);
    REQUIRE(std::abs(outd[7] - 4852140.073) < 1e-3);
    REQUIRE(std::abs(outt[0] - 1331.6311) < 1e-3);
    REQUIRE(std::abs(outt[1] - 1321.5417) < 1e-3);
  }

  f107 = 150.;
  ap[0] = 40;
  ap[1] = 0;
  SUBCASE("Example 10, testing for precesion 1e-3.") {
    msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
           outt);
    REQUIRE(std::abs(outd[0] - 557795.943) < 1e-3);
    REQUIRE(std::abs(outd[1] - 111566132.416) < 1e-3);
    REQUIRE(std::abs(outd[2] - 38317548.571) < 1e-3);
    REQUIRE(std::abs(outd[3] - 1717539.631) < 1e-3);
    REQUIRE(std::abs(outd[4] - 12488.943) < 1e-3);
    REQUIRE(std::abs(outd[5] - 4.91e-15) < 1e-17);
    REQUIRE(std::abs(outd[6] - 26782.682) < 1e-3);
    REQUIRE(std::abs(outd[7] - 3176660.134) < 1e-3);
    REQUIRE(std::abs(outt[0] - 1384.0740) < 1e-3);
    REQUIRE(std::abs(outt[1] - 1373.4963) < 1e-3);
  }
}