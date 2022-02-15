#include "datetime/datetime_write.hpp"
#include "datetime/dtcalendar.hpp"
#include "planetpos.hpp"
#include <cassert>
#include <cstdio>
#include <random>

int num_tests = 1000;
std::random_device rd;  // obtain a random number from hardware
std::mt19937 gen(rd()); // seed the generator
std::uniform_int_distribution<> ydistr(1995, 2035); // define the range
std::uniform_int_distribution<> mdistr(1, 12);       // define the range
std::uniform_int_distribution<> ddistr(1, 28);       // define the range

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s [DE/SPK KERNEL] [LEAPSEC/LSK KERNEL]\n",
            argv[0]);
    return 1;
  }
  // to compute the planet's position via cspice, we need to load:
  // 1. the planetary ephemeris (SPK) kernel
  // 2. the leap-second (aka LSK) kernel
  dso::cspice::load_if_unloaded_spk(argv[1]);
  dso::cspice::load_if_unloaded_lsk(argv[2]);

  double vmoon1[3], vmoon2[3], vmoon_ref[3];
  double mr1=0e0, mr2=0e0;
  double rdiff;

  for (int i = 0; i < num_tests; i++) {
    dso::datetime<dso::seconds> t(
        dso::year(ydistr(gen)), dso::month(mdistr(gen)),
        dso::day_of_month(ddistr(gen)), dso::seconds(0));

    dso::moon_vector_cspice(t, vmoon_ref);
    printf("CSPICE/ANIF: %+15.1f %+15.1f %+15.1f\n", vmoon_ref[0], vmoon_ref[1],
           vmoon_ref[2]);

    dso::moon_vector_approx(t, vmoon1);
    printf("Monetnbruck: %+15.1f %+15.1f %+15.1f %15.1f %15.1f %15.1f\n",
           vmoon1[0], vmoon1[1], vmoon1[2], std::abs(vmoon1[0] - vmoon_ref[0]),
           std::abs(vmoon1[1] - vmoon_ref[1]), std::abs(vmoon1[2] - vmoon_ref[2]));
    rdiff = std::sqrt((vmoon1[0] - vmoon_ref[0]) * (vmoon1[0] - vmoon_ref[0]) +
                      (vmoon1[1] - vmoon_ref[1]) * (vmoon1[1] - vmoon_ref[1]) +
                      (vmoon1[2] - vmoon_ref[2]) * (vmoon1[2] - vmoon_ref[2]));
    mr1 += (1/(i+1)) * (rdiff - mr1);

    dso::moon_vector_vallado(t, vmoon2);
    printf("Vallado    : %+15.1f %+15.1f %+15.1f %15.1f %15.1f %15.1f\n",
           vmoon2[0], vmoon2[1], vmoon2[2], std::abs(vmoon2[0] - vmoon_ref[0]),
           std::abs(vmoon2[1] - vmoon_ref[1]), std::abs(vmoon2[2] - vmoon_ref[2]));
    rdiff = std::sqrt((vmoon2[0] - vmoon_ref[0]) * (vmoon2[0] - vmoon_ref[0]) +
                      (vmoon2[1] - vmoon_ref[1]) * (vmoon2[1] - vmoon_ref[1]) +
                      (vmoon2[2] - vmoon_ref[2]) * (vmoon2[2] - vmoon_ref[2]));
    mr2 += (1/(i+1)) * (rdiff - mr2);

    //dso::moon_vector_curtis(t, vmoon4);
    //printf("Curtis     : %+15.1f %+15.1f %+15.1f %15.1f %15.1f %15.1f\n",
    //       vmoon4[0], vmoon4[1], vmoon4[2], std::abs(vmoon4[0] - vmoon_ref[0]),
    //       std::abs(vmoon4[1] - vmoon_ref[1]), std::abs(vmoon4[2] - vmoon_ref[2]));
    //rdiff = std::sqrt((vmoon4[0] - vmoon_ref[0]) * (vmoon4[0] - vmoon_ref[0]) +
    //                  (vmoon4[1] - vmoon_ref[1]) * (vmoon4[1] - vmoon_ref[1]) +
    //                  (vmoon4[2] - vmoon_ref[2]) * (vmoon4[2] - vmoon_ref[2]));
    //mr3 += (1/(i+1)) * (rdiff - mr3);
  }

  printf("Distance Average Differences (km)\n");
  printf("Montenbruck: %10.1f\n", mr1);
  printf("Vallado    : %10.1f\n", mr2);
  //printf("Curtis     : %10.1f\n", mr3);

  dso::datetime<dso::seconds> tb(
        dso::year(1994), dso::month(4),
        dso::day_of_month(28), dso::seconds(0));
  double targ = - 0.056796714;
  dso::moon_vector_vallado(targ, vmoon2);
  dso::moon_vector_cspice(tb, vmoon_ref);
  printf("Note, book example: %15.1f %15.1f %15.1f\n", vmoon2[0],vmoon2[1],vmoon2[2]);
  printf("CSPICE            : %15.1f %15.1f %15.1f\n", vmoon_ref[0],vmoon_ref[1],vmoon_ref[2]);
  
  dso::datetime<dso::seconds> tc(
        dso::year(2013), dso::month(7),
        dso::day_of_month(25), dso::seconds(8*3600L));
  dso::moon_vector_vallado(tc, vmoon2);
  dso::moon_vector_cspice(tc, vmoon_ref);
  printf("Note, book example: %15.1f %15.1f %15.1f\n", vmoon2[0],vmoon2[1],vmoon2[2]);
  printf("CSPICE            : %15.1f %15.1f %15.1f\n", vmoon_ref[0],vmoon_ref[1],vmoon_ref[2]);

  return 0;
}
