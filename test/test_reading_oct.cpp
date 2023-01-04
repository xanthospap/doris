#include "tides.hpp"
#include <cstdio>

constexpr const int Degree = 5;
constexpr const int Order = Degree;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <OCEAN TIDE FILE>\n", argv[0]);
    return 1;
  }

  std::vector<dso::DoodsonOceanTideConstituent> vdds;
  // if (dso::inspect_octide_coefficients(argv[1], vdds)) {
  //   fprintf(stderr, "Failed reading inpit file!\n");
  //   return 1;
  // }
  if (dso::memmap_octide_coefficients(argv[1], vdds, Degree, Order, 3)) {
    fprintf(stderr, "Failed reading inpit file!\n");
    return 1;
  }

  int i = 0;
  for (const auto &d : vdds) {
    char buf[32];
    printf("#%02d Doodson: %s -> maxl:%d maxm:%d\n", i,
           d.doodson_number().str(buf), d.max_degree(), d.max_order());
    ++i;
  }

  return 0;
}
