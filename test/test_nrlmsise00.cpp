#include "atmosphere.hpp"
#include <cstdio>

using dso::air_density_models::nrlmsise00::nrlmsise_flags;
using dso::air_density_models::nrlmsise00::nrlmsise_input;
using dso::air_density_models::nrlmsise00::nrlmsise_output;
using dso::air_density_models::nrlmsise00::gtd7;

void test_gtd7() {
  nrlmsise_output output[17];
  nrlmsise_input input[17];
  nrlmsise_flags flags;

  /* input values */
  //for (int i = 0; i < 7; i++)
  //  aph.a[i] = 100;
  
  flags.switches[0] = 0;
  for (int i = 1; i < 24; i++)
    flags.switches[i] = 1;
  
  for (int i = 0; i < 17; i++) {
    input[i].doy = 172;
    input[i].year = 0; /* without effect */
    input[i].sec = 29000;
    input[i].alt = 400;
    input[i].g_lat = 60;
    input[i].g_lon = -70;
    input[i].lst = 16;
    input[i].f107A = 150;
    input[i].f107 = 150;
    input[i].ap = 4;
    //for (int j = 0; j < 7; j++)
    //  input[i].ap_array[j] = 100e0;
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
  input[4].g_lat = 0;
  input[5].g_lon = 0;
  input[6].lst = 4;
  input[7].f107A = 70;
  input[8].f107 = 180;
  input[9].ap = 40;
  for (int i=15; i<17; i++) for (int j=0; j<7; j++) input[i].ap_array[j] = 100e0;
  // input[15].ap_a = aph;
  // input[16].ap_a = aph;
  
  /* evaluate 0 to 14 */
  for (int i = 0; i < 15; i++)
    gtd7(input[i], flags, output[i]);
  
  /* evaluate 15 and 16 */
  flags.switches[9] = -1;
  for (int i = 15; i < 17; i++)
    gtd7(input[i], flags, output[i]);
  
  /* output type 1 */
  for (int i = 0; i < 17; i++) {
    printf("\n");
    for (int j = 0; j < 9; j++)
      printf("%E ", output[i].d[j]);
    printf("%E ", output[i].t[0]);
    printf("%E \n", output[i].t[1]);
    /* DL omitted */
  }

  /* output type 2 */
  int i,j;
  for (i = 0; i < 3; i++) {
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
      printf("         %3.0f", input[i * 5 + j].g_lat);
    printf("\nLONG  ");
    for (j = 0; j < 5; j++)
      printf("         %3.0f", input[i * 5 + j].g_lon);
    printf("\nLST   ");
    for (j = 0; j < 5; j++)
      printf("       %5.0f", input[i * 5 + j].lst);
    printf("\nF107A ");
    for (j = 0; j < 5; j++)
      printf("         %3.0f", input[i * 5 + j].f107A);
    printf("\nF107  ");
    for (j = 0; j < 5; j++)
      printf("         %3.0f", input[i * 5 + j].f107);
    printf("\n\n");
    printf("\nTINF  ");
    for (j = 0; j < 5; j++)
      printf("     %7.2f", output[i * 5 + j].t[0]);
    printf("\nTG    ");
    for (j = 0; j < 5; j++)
      printf("     %7.2f", output[i * 5 + j].t[1]);
    printf("\nHE    ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", output[i * 5 + j].d[0]);
    printf("\nO     ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", output[i * 5 + j].d[1]);
    printf("\nN2    ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", output[i * 5 + j].d[2]);
    printf("\nO2    ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", output[i * 5 + j].d[3]);
    printf("\nAR    ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", output[i * 5 + j].d[4]);
    printf("\nH     ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", output[i * 5 + j].d[6]);
    printf("\nN     ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", output[i * 5 + j].d[7]);
    printf("\nANM 0 ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", output[i * 5 + j].d[8]);
    printf("\nRHO   ");
    for (j = 0; j < 5; j++)
      printf("   %1.3e", output[i * 5 + j].d[5]);
    printf("\n");
  }
  printf("\n");
}

int main() {
  test_gtd7();
  return 0;
}