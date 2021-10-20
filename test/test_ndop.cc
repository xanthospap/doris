#include "doris_rinex.hpp"
#include <iostream>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "\nUsage: %s [RINEX_FILE]\n", argv[0]);
    return 1;
  }

  ids::DorisObsRinex rnx(argv[1]);
  int error = rnx.get_doppler_counts();

  printf("Extract Doppler counts returned %d\n", error);
  return 0;
}
