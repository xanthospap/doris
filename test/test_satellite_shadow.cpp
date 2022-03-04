#include "satellite.hpp"
#include "planetpos.hpp"
#include "sp3/sp3.hpp"
#include <cstdio>

void shadow_coeff(const dso::Sp3DataBlock *block, const double *vsun) noexcept {
  dso::Vector3 rsat{
      {block->state[0] * 1e3, block->state[1] * 1e3, block->state[2] * 1e3}};
  dso::Vector3 rsun{{vsun[0] * 1e3, vsun[1] * 1e3, vsun[2] * 1e3}};
  double c1 = utest::montebruck_shadow(rsat, rsun);
  double c2 = utest::conical_shadow(rsat, rsun);
  printf("%.5f f1=%.3f f2=%.3f\n", block->t.as_mjd(), c1, c2);
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
    fprintf(stderr, "Usage: %s [SP3] [DE/SPK KERNEL] [LEAPSEC/LSK KERNEL]\n",
            argv[0]);
    return 1;
  }

  dso::Sp3c sp3(argv[1]);
  #ifdef DEBUG
  sp3.print_members();
  #endif
  if (sp3.num_sats() != 1) {
      fprintf(stderr, "More than one satellites found in SP3 file!\n");
      return 1;
  }
  dso::sp3::SatelliteId sv;
  sv.set_id(sp3.sattellite_vector()[0].id);
  dso::Sp3DataBlock block;

  // to compute the planet's position via cspice, we need to load:
  // 1. the planetary ephemeris (SPK) kernel
  // 2. the leap-second (aka LSK) kernel
  dso::cspice::load_if_unloaded_spk(argv[2]);
  dso::cspice::load_if_unloaded_lsk(argv[3]);

    double vsun[3];
     // let's try reading the records; note that -1 denotes EOF
  int j;
  std::size_t rec_count = 0;
  do {
    j = sp3.get_next_data_block(sv, block);
    if (j > 0) {
      fprintf(stderr, "Something went wrong ....status = %3d\n", j);
      return 1;
    } else if (j == -1) {
      printf("EOF encountered; Sp3 file read through!\n");
    }

    // sun position vector, geocentric [km]
    dso::sun_vector_cspice(block.t, vsun);
    
    // compute shadow functions ...
    shadow_coeff(&block, vsun);

    ++rec_count;
  } while (!j);

  printf("Num of records read: %6lu\n", rec_count);

}