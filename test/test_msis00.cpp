#include "atmosphere.hpp"
#include <cstdio>

using dso::air_density_models::nrlmsise00::gtd7;

struct auxin {
    double doy;
    int year;
    double sec;
    double alt;
    double lat;
    double lon;
    double lst;
    double f107;
    double f107A;
    double ap;
};

struct auxout {
    double outd[10];
    double outt[2];
};

int main() {
  // input values
  double apha[7];
  for (int i = 0; i < 7; i++) apha[i] = 100;
  
    int switches[24];
  switches[0] = 0;
  for (int i = 1; i < 24; i++)
    switches[i] = 1;
  
  auxin input[17];
  for (int i = 0; i < 17; i++) {
    input[i].doy = 172;
    input[i].year = 0; // without effect
    input[i].sec = 29000;
    input[i].alt = 400;
    input[i].lat = 60;
    input[i].lon = -70;
    input[i].lst = 16;
    input[i].f107A = 150;
    input[i].f107 = 150;
    input[i].ap = 4;
  }

  input[1].doy = 81;
  input[2].sec = 75000;
  input[2].alt = 1000;
  input[3].alt = 100;
  input[10].alt = 0;
  input[11].alt = 10;
  input[12].alt = 30;
  input[13].alt = 50;
  input[14].alt = 70;
  input[16].alt = 100;
  input[4].lat = 0;
  input[5].lon = 0;
  input[6].lst = 4;
  input[7].f107A = 70;
  input[8].f107 = 180;
  input[9].ap = 40;
  // input[15].ap_a = &aph;
  // input[16].ap_a = &aph;

  // evaluate 0 to 14
  auxout out[17];
  for (int i = 0; i < 15; i++)
    gtd7(switches, input[i].doy, input[i].sec,
         input[i].lat, input[i].lon, input[i].lst, input[i].f107, input[i].f107A,
         input[i].alt, apha, out[i].outd, out[i].outt);

  // evaluate 15 and 16
  switches[9] = -1;
  for (int i = 15; i < 17; i++)
    gtd7(switches, input[i].doy, input[i].sec,
         input[i].lat, input[i].lon, input[i].lst, input[i].f107, input[i].f107A,
         input[i].alt, apha, out[i].outd, out[i].outt);
  
  // output type 1
  for (int i = 0; i < 17; i++) {
    printf("\n");
    for (int j = 0; j < 9; j++)
      printf("%E ", out[i].outd[j]);
    
    printf("%E ", out[i].outt[0]);
    printf("%E \n", out[i].outt[1]);
  }
}