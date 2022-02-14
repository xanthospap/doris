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

  double vsun1[3], vsun2[3], vsun_ref[3], vsun4[3];
  double mr1=0e0, mr2=0e0, mr3=0e0;
  double rdiff;

  for (int i = 0; i < num_tests; i++) {
    dso::datetime<dso::seconds> t(
        dso::year(ydistr(gen)), dso::month(mdistr(gen)),
        dso::day_of_month(ddistr(gen)), dso::seconds(0));

    dso::sun_vector_cspice(t, vsun_ref);
    printf("CSPICE/ANIF: %+15.1f %+15.1f %+15.1f\n", vsun_ref[0], vsun_ref[1],
           vsun_ref[2]);

    dso::sun_vector_montenbruck(t, vsun1);
    printf("Monetnbruck: %+15.1f %+15.1f %+15.1f %15.1f %15.1f %15.1f\n",
           vsun1[0], vsun1[1], vsun1[2], std::abs(vsun1[0] - vsun_ref[0]),
           std::abs(vsun1[1] - vsun_ref[1]), std::abs(vsun1[2] - vsun_ref[2]));
    rdiff = std::sqrt((vsun1[0] - vsun_ref[0]) * (vsun1[0] - vsun_ref[0]) +
                      (vsun1[1] - vsun_ref[1]) * (vsun1[1] - vsun_ref[1]) +
                      (vsun1[2] - vsun_ref[2]) * (vsun1[2] - vsun_ref[2]));
    mr1 += (1/(i+1)) * (rdiff - mr1);

    dso::sun_vector_vallado(t, vsun2);
    printf("Vallado    : %+15.1f %+15.1f %+15.1f %15.1f %15.1f %15.1f\n",
           vsun2[0], vsun2[1], vsun2[2], std::abs(vsun2[0] - vsun_ref[0]),
           std::abs(vsun2[1] - vsun_ref[1]), std::abs(vsun2[2] - vsun_ref[2]));
    rdiff = std::sqrt((vsun2[0] - vsun_ref[0]) * (vsun2[0] - vsun_ref[0]) +
                      (vsun2[1] - vsun_ref[1]) * (vsun2[1] - vsun_ref[1]) +
                      (vsun2[2] - vsun_ref[2]) * (vsun2[2] - vsun_ref[2]));
    mr2 += (1/(i+1)) * (rdiff - mr2);

    dso::sun_vector_curtis(t, vsun4);
    printf("Curtis     : %+15.1f %+15.1f %+15.1f %15.1f %15.1f %15.1f\n",
           vsun4[0], vsun4[1], vsun4[2], std::abs(vsun4[0] - vsun_ref[0]),
           std::abs(vsun4[1] - vsun_ref[1]), std::abs(vsun4[2] - vsun_ref[2]));
    rdiff = std::sqrt((vsun4[0] - vsun_ref[0]) * (vsun4[0] - vsun_ref[0]) +
                      (vsun4[1] - vsun_ref[1]) * (vsun4[1] - vsun_ref[1]) +
                      (vsun4[2] - vsun_ref[2]) * (vsun4[2] - vsun_ref[2]));
    mr3 += (1/(i+1)) * (rdiff - mr3);
  }

  printf("Distance Average Differences (km)\n");
  printf("Montenbruck: %10.1f\n", mr1);
  printf("Vallado    : %10.1f\n", mr2);
  printf("Curtis     : %10.1f\n", mr3);

  return 0;
}
