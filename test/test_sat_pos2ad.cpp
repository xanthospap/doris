#include "satellite.hpp"
#include "geodesy/units.hpp"
#include <cstdio>

// Example: calculate right ascension and declination from position vector
int main() {
  // Example 4.1 for Curtis
  dso::Vector3 r {-5368e0, -1784e0, 3691e0}; // position vector in [km]

  double ra, dec;
  dso::pos2ad(r, ra, dec);

  printf("Declination    : %+10.3f\n", dso::rad2deg(dec));
  printf("Right Ascension: %+10.3f\n", dso::rad2deg(ra));

  // results reported from the book:
  // δ = 33.12 [deg]
  // α = 198.4 [deg]
  printf("Delta-Declination : %5.2f\n", dso::rad2deg(dec)-33.12e0);
  printf("Delta-R. Ascension: %5.2f\n", dso::rad2deg(ra)-198.4e0);

  return 0;
}
