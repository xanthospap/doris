#include "tides.hpp"
#include <cstdio>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <OCEAN TIDE FILE>\n", argv[0]);
    return 1;
  }

  std::vector<dso::DoodsonOceanTideConstituent> vdds;
  if (dso::inspect_octide_coefficients(argv[1], vdds)) {
    fprintf(stderr, "Failed reading inpit file!\n");
    return 1;
  }

  int i = 0;
  for (const auto &d : vdds) {
    char buf[32];
    sprintf(buf, "%d%+2d%+2d.%+2d%+2d%+2d", d.doodson[0], d.doodson[1], d.doodson[2],
            d.doodson[3], d.doodson[4], d.doodson[5]);
    printf("#%02d Doodson: %s minl:%d maxl:%d maxm:%d\n", i, buf, d.minl, d.maxl,
           d.maxm);
    ++i;
  }

  return 0;
}
