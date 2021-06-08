#include "antenna_pcv.hpp"
#include <cassert>
#include <cstdio>

using namespace ids;

int main() {
  static_assert(AntennaOffset<GroundAntennaType::Starec_B, 1>::dnorth() == 0e0);

  static_assert(AntennaOffset<GroundAntennaType::Alcatel, 1>::dnorth() == 0e0);

  static_assert(AntennaOffset<GroundAntennaType::Alcatel, 1>::dup() == 510.0e0);

  static_assert(AntennaOffset<GroundAntennaType::Starec_B, 1>::dnorth() ==
                AntennaOffset<GroundAntennaType::Starec_C, 1>::dnorth());

  int d1_out_of_bounds, d2_out_of_bounds;
  printf("Zenith      Alcatel_1  Alcatel_2  Starec_1  Starec_2\n");
  for (double zenith = 0e0; zenith <= 90e0; zenith += 0.3e0) {
    double zrad = ngpt::deg2rad<double>(zenith);
    double alc1 = AntennaOffset<GroundAntennaType::Alcatel, 1>::pcv(
        zrad, d1_out_of_bounds);
    double alc2 = AntennaOffset<GroundAntennaType::Alcatel, 2>::pcv(
        zrad, d2_out_of_bounds);
    assert(d1_out_of_bounds + d2_out_of_bounds == 0);
    /*
    if (d1_out_of_bounds || d2_out_of_bounds) {
      fprintf(stderr, "[ERROR] Out-of-bounds zenith angle!");
      return 2;
    }*/
    double str1 = AntennaOffset<GroundAntennaType::Starec_B, 1>::pcv(
        zrad, d1_out_of_bounds);
    double str2 = AntennaOffset<GroundAntennaType::Starec_B, 2>::pcv(
        zrad, d2_out_of_bounds);
    assert(d1_out_of_bounds + d2_out_of_bounds == 0);
    /*
    if (d1_out_of_bounds || d2_out_of_bounds) {
      fprintf(stderr, "[ERROR] Out-of-bounds zenith angle!");
      return 2;
    }*/
    printf("%+10.4f %+10.4f %+10.4f %+10.4f %+10.4f\n", zenith, alc1, alc2,
           str1, str2);
  }

  return 0;
}
