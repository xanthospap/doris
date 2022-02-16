#include "geodesy/units.hpp"
#include "satellite.hpp"
#include <cstdio>

// Example: calculate orbital elements from state vector
int main() {
  // Example 4.3 from Curtis
  double state[] = {-6045e3,  -3490e3, 2500e3,
                    -3.457e3, 6.618e3, 2.533e3}; // [m] and [m/s]
  double elements[6];
  const double answer[] = {58310e0,  153.2e0, 255.3e0,
                           0.1712e0, 20.07e0, 28.45e0}; // [km]

  dso::state2kepler(state, elements);
  
  // for comparisson, transform radians to degrees ... and km^2 to m^2
  double h = elements[0] * 1e-6;
  double i = dso::rad2deg(elements[1]);
  double Omega = dso::rad2deg(elements[2]);
  double omega = dso::rad2deg(elements[4]);
  double theta = dso::rad2deg(elements[5]);
  // and print results
  printf("State Vector to Keplerian elements results:\n");
  printf("h    : %10.3f Diff %10.1f\n", h, h - answer[0]);
  printf("i    : %10.3f Diff %10.2f\n", i, i - answer[1]);
  printf("Omega: %10.3f Diff %10.2f\n", Omega, Omega - answer[2]);
  printf("e    : %10.3f Diff %10.4f\n", elements[3], elements[3] - answer[3]);
  printf("omega: %10.3f Diff %10.2f\n", omega, omega - answer[4]);
  printf("theta: %10.3f Diff %10.2f\n", theta, theta - answer[5]);

  // if we wanted to perform the computations in [km] (instead of [m]) ...
  double state3[6], elements3[6];
  std::memcpy(state3, state, sizeof(double)*6);
  for (int j = 0; j < 6; j++)
    state3[j] *= 1e3; // [m] to [km]
  dso::state2kepler(state3, elements3, iers2010::GMe / 1e9);
  // for comparisson, transform radians to degrees
  elements3[1] = dso::rad2deg(elements3[1]);
  elements3[2] = dso::rad2deg(elements3[2]);
  elements3[4] = dso::rad2deg(elements3[4]);
  elements3[5] = dso::rad2deg(elements3[5]);
  // print results
  printf("h    : %10.3f Diff %10.1f\n", elements3[0], elements3[0] - answer[0]);
  printf("i    : %10.3f Diff %10.2f\n", elements3[1], elements3[1] - answer[1]);
  printf("Omega: %10.3f Diff %10.2f\n", elements3[2], elements3[2] - answer[2]);
  printf("e    : %10.3f Diff %10.4f\n", elements3[3], elements3[3] - answer[3]);
  printf("omega: %10.3f Diff %10.2f\n", elements3[4], elements3[4] - answer[4]);
  printf("theta: %10.3f Diff %10.2f\n", elements3[5], elements3[5] - answer[5]);

  // let's now use the inverse computation to go back, aka from orbital 
  // elements to state vector
  double state2[6];
  dso::kepler2state(elements, state2);
  printf("Keplerian elements to State vector results:\n");
  printf("X : %15.1f (diff from input: %+15.12f)\n", state2[0], state2[0]-state[0]);
  printf("Y : %15.1f (diff from input: %+15.12f)\n", state2[1], state2[1]-state[1]);
  printf("Z : %15.1f (diff from input: %+15.12f)\n", state2[2], state2[2]-state[2]);
  printf("DX: %15.1f (diff from input: %+15.12f)\n", state2[3], state2[3]-state[3]);
  printf("DY: %15.1f (diff from input: %+15.12f)\n", state2[4], state2[4]-state[4]);
  printf("SZ: %15.1f (diff from input: %+15.12f)\n", state2[5], state2[5]-state[5]);

  return 0;
}
