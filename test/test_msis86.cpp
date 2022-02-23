#include "atmosphere.hpp"
#include <cstdio>
#include <cstring>

using dso::air_density_models::msis86::msis86;

int main() {
  int day = 172;       // day of the year
  double sec = 29000.; // utc in seconds
  double alt = 400.;   // altitude in km
  double glat = 60.;   // latitude geodetic in degrees
  double glong = -70.; // longitude in degrees
  double stl = 16.;    // local solar time in hours
  double f107a = 150.; // mean solar flux over 90 days
  double f107 = 150.;  // solar flux from previous day
  double ap[] = {4, 4, 4, 4, 4, 4, 4};
  int switches[24];
  std::memset(switches, 1, sizeof(int) * 24);
  /*
  dayly Ap
  3 hr ap index for current time
  3 hr ap index for 3 hrs before current time
  3 hr ap index for 6 hrs before current time
  3 hr ap index for 9 hrs before current time
  average of eight 3 hr ap indicies from 12 to 33 hrs prior to current time
  average of eight 3 hr ap indicies from 36 to 59 hrs prior to current time
  */

  double outd[8], outt[2];
  msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
         outt);
  for (int i = 0; i < 8; i++)
    printf("%.2E ", outd[i]);
  printf("\n    %.3E %.3E \n", outt[0], outt[1]);

  //day = 81;
  //msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
  //       outt);
  //for (int i = 0; i < 8; i++)
  //  printf("%.2E ", outd[i]);
  //printf("\n    %.3E %.3E \n", outt[0], outt[1]);

  //day = 172;
  //sec = 75000.;
  //msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
  //       outt);
  //for (int i = 0; i < 8; i++)
  //  printf("%.2E ", outd[i]);
  //printf("\n    %.3E %.3E \n", outt[0], outt[1]);

  //sec = 29000.;
  //alt = 100.;
  //ap[6] = 40;
  //msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
  //       outt);
  //for (int i = 0; i < 8; i++)
  //  printf("%.2E ", outd[i]);
  //printf("\n    %.3E %.3E \n", outt[0], outt[1]);

  //alt = 400.;
  //glat = 0.;
  //ap[5] = 40;
  //ap[6] = 0;
  //msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
  //       outt);
  //for (int i = 0; i < 8; i++)
  //  printf("%.2E ", outd[i]);
  //printf("\n    %.3E %.3E \n", outt[0], outt[1]);

  //glat = 60.;
  //glong = 0.;
  //ap[4] = 40;
  //ap[5] = 0;
  //msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
  //       outt);
  //for (int i = 0; i < 8; i++)
  //  printf("%.2E ", outd[i]);
  //printf("\n    %.3E %.3E \n", outt[0], outt[1]);

  //glong = -70.;
  //stl = 4.;
  //ap[3] = 40;
  //ap[4] = 0;
  //msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
  //       outt);
  //for (int i = 0; i < 8; i++)
  //  printf("%.2E ", outd[i]);
  //printf("\n    %.3E %.3E \n", outt[0], outt[1]);

  //stl = 16.;
  //f107a = 70.;
  //ap[2] = 40;
  //ap[3] = 0;
  //msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
  //       outt);
  //for (int i = 0; i < 8; i++)
  //  printf("%.2E ", outd[i]);
  //printf("\n    %.3E %.3E \n", outt[0], outt[1]);

  //f107a = 150.;
  //f107 = 180.;
  //ap[1] = 40;
  //ap[2] = 0;
  //msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
  //       outt);
  //for (int i = 0; i < 8; i++)
  //  printf("%.2E ", outd[i]);
  //printf("\n    %.3E %.3E \n", outt[0], outt[1]);

  //f107 = 150.;
  //ap[0] = 40;
  //ap[1] = 0;
  //msis86(day, sec, alt, glat, glong, stl, f107a, f107, switches, ap, outd,
  //       outt);
  //for (int i = 0; i < 8; i++)
  //  printf("%.2E ", outd[i]);
  //printf("\n    %.3E %.3E \n", outt[0], outt[1]);

  return 0;
}