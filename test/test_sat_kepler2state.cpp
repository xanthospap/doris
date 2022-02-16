#include "geodesy/units.hpp"
#include "satellite.hpp"
#include <cstdio>

using dso::deg2rad;
using dso::rad2deg;

// Example: calculate orbital elements from state vector
int main() {
  // Example 4.7 from Curtis
  const double orbelements[] = {80'000e0,      deg2rad(30e0),
                                deg2rad(40e0), 1.4e0,
                                deg2rad(60e0), deg2rad(30e0)}; // [km] and [rad]
  const double answer[] = {-4'040e0, 4'815e0, 3'629e0,
                           -10.39e0, -4.772e0, 1.744e0}; // [km] and [km/s]
  double state[6];

  dso::kepler2state(orbelements, state,
                    iers2010::GMe / 1e9); // perform calculations in [km]
  printf("Keplerian elements to State vector results:\n");
  printf("X : %15.1f (diff from input: %+15.3f)\n", state[0],
         state[0] - answer[0]);
  printf("Y : %15.1f (diff from input: %+15.3f)\n", state[1],
         state[1] - answer[1]);
  printf("Z : %15.1f (diff from input: %+15.3f)\n", state[2],
         state[2] - answer[2]);
  printf("DX: %15.4f (diff from input: %+15.4f)\n", state[3],
         state[3] - answer[3]);
  printf("DY: %15.4f (diff from input: %+15.4f)\n", state[4],
         state[4] - answer[4]);
  printf("SZ: %15.4f (diff from input: %+15.4f)\n", state[5],
         state[5] - answer[5]);

  double ele[6];
  dso::state2kepler(state, ele,
                    iers2010::GMe / 1e9); // perform calculations in [km]
  // and print results
  printf("State Vector to Keplerian elements results:\n");
  printf("h    : %10.3f Diff %10.1f\n", ele[0], ele[0] - orbelements[0]);
  printf("i    : %10.3f Diff %10.2f\n", rad2deg(ele[1]),
         rad2deg(ele[1] - orbelements[1]));
  printf("Omega: %10.3f Diff %10.2f\n", rad2deg(ele[2]),
         rad2deg(ele[2] - orbelements[2]));
  printf("e    : %10.3f Diff %10.4f\n", ele[3], ele[3] - orbelements[3]);
  printf("omega: %10.3f Diff %10.2f\n", rad2deg(ele[4]),
         rad2deg(ele[4] - orbelements[4]));
  printf("theta: %10.3f Diff %10.2f\n", rad2deg(ele[5]),
         rad2deg(ele[5] - orbelements[5]));

  return 0;
}
