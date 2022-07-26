#include "atmosphere.hpp"
#include <cassert>
#include <cstdio>

#define DPRECISION 1e-12
#define TPRECISION 1e-9
#define TOTALD_PRECISION 1e-15

int main() {

  dso::nrlmsise00::ApArray aph;
  for (int i = 0; i < 7; i++)
    aph.a[i] = 100;

  dso::nrlmsise00::InParams input[17];
  dso::nrlmsise00::OutParams out[17];

  for (int i = 0; i < 17; i++) {
    input[i].doy = 172;     // day of year
    input[i].year = 0;      // without effect
    input[i].sec = 29000e0; // ut
    input[i].alt = 400e0;
    input[i].glat = 60e0;
    input[i].glon = -70e0;
    input[i].lst = 16e0;
    input[i].f107A = 150e0;
    input[i].f107 = 150e0;
    input[i].ap = 4e0; // magnetic_index
    input[i].aparr = aph;
  }

  input[1].doy = 81;
  input[2].sec = 75000e0;
  input[2].alt = 1000e0;
  input[3].alt = 100e0;
  input[10].alt = 0e0;
  input[11].alt = 10e0;
  input[12].alt = 30e0;
  input[13].alt = 50e0;
  input[14].alt = 70e0;
  input[16].alt = 100e0;
  input[4].glat = 0e0;
  input[5].glon = 0e0;
  input[6].lst = 4e0;
  input[7].f107A = 70e0;
  input[8].f107 = 180e0;
  input[9].ap = 40e0;
  // input[15].ap_a = &aph;
  // input[16].ap_a = &aph;

  dso::Nrlmsise00 Msise;
  int mass = 48;

  // evaluate 0 to 14
  for (int i = 0; i < 15; i++) {
    input[i].set_switches_on();
    Msise.gtd7(&input[i], mass, &out[i]);
    return 0;
  }

  // evaluate 15 and 16
  for (int i = 15; i < 17; i++) {
    Msise.gtd7(&input[i], mass, &out[i]);
  }

  /* output type 2 */
  int j;
  for (int i = 0; i < 3; i++) {
    printf("\n");
    printf("\nDAY   ");
    for (j = 0; j < 5; j++)
      printf("         %3i", input[i * 5 + j].doy);
    printf("\nUT    ");
    for (j = 0; j < 5; j++)
      printf("       %5.0f", input[i * 5 + j].sec);
    printf("\nALT   ");
    for (j = 0; j < 5; j++)
      printf("        %4.0f", input[i * 5 + j].alt);
    printf("\nLAT   ");
    for (j = 0; j < 5; j++)
      printf("         %3.0f", input[i * 5 + j].glat);
    printf("\nLONG  ");
    for (j = 0; j < 5; j++)
      printf("         %3.0f", input[i * 5 + j].glon);
    printf("\nLST   ");
    for (j = 0; j < 5; j++)
      printf("       %5.0f", input[i * 5 + j].lst);
    printf("\nF107A ");
    for (j = 0; j < 5; j++)
      printf("         %3.0f", input[i * 5 + j].f107A);
    printf("\nF107  ");
    for (j = 0; j < 5; j++)
      printf("         %3.0f", input[i * 5 + j].f107);
    printf("\nAP    ");
    for (j = 0; j < 5; j++)
      printf("         %3.0f", input[i * 5 + j].ap);
    printf("\n\n");
    printf("\nTINF  ");
    for (j = 0; j < 5; j++)
      printf("     %7.2f", out[i * 5 + j].t[0]);
    printf("\nTG    ");
    for (j = 0; j < 5; j++)
      printf("     %7.2f", out[i * 5 + j].t[1]);
    printf("\nHE    ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", out[i * 5 + j].d[0]);
    printf("\nO     ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", out[i * 5 + j].d[1]);
    printf("\nN2    ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", out[i * 5 + j].d[2]);
    printf("\nO2    ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", out[i * 5 + j].d[3]);
    printf("\nAR    ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", out[i * 5 + j].d[4]);
    printf("\nH     ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", out[i * 5 + j].d[6]);
    printf("\nN     ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", out[i * 5 + j].d[7]);
    printf("\nANM 0 ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", out[i * 5 + j].d[8]);
    printf("\nRHO   ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", out[i * 5 + j].d[5]);
    printf("\n");
  }
  printf("\n");
}
