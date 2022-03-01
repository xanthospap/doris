#include "atmosphere.hpp"
#include <cassert>
#include <cmath>
#ifdef DEBUG
#include <cstdio>
#endif

struct AtmLayer {
  int limit_low;          // [km]
  int limit_high;         // [km]
  int base_altitude;      // [km]
  double nominal_density; // [kg/m3]
  double scale_height;    // [km]
};

const AtmLayer AtmTable[] = {
    {0, 25, 0, 1.225e0, 7.249e0},
    {25, 30, 25, 3.899e-2, 6.349e0},
    {30, 40, 30, 1.774e-2, 6.682e0},
    {40, 50, 40, 3.972e-3, 7.554e0},
    {50, 60, 50, 1.057e-3, 8.382e0},
    {60, 70, 60, 3.206e-4, 7.714e0},
    {70, 80, 70, 8.770e-5, 6.549e0},
    {80, 90, 80, 1.905e-5, 5.799e0},
    {90, 100, 90, 3.396e-6, 5.382e0},
    {100, 110, 100, 5.297e-7, 5.877e0},
    {110, 120, 110, 9.661e-8, 7.263e0},
    {120, 130, 120, 2.438e-8, 9.473e0},
    {130, 140, 130, 8.484e-9, 12.636e0},
    {140, 150, 140, 3.845e-9, 16.149e0},
    {150, 180, 150, 2.070e-9, 22.523e0},
    {180, 200, 180, 5.464e-10, 29740e0},
    {200, 250, 200, 2.789e-10, 37.105e0},
    {250, 300, 250, 7.248e-11, 45.546e0},
    {300, 350, 300, 2.418e-11, 53.628e0},
    {350, 400, 350, 9.518e-12, 53.298e0},
    {400, 450, 400, 3.725e-12, 58.515e0},
    {450, 500, 450, 1.585e-12, 60.828e0},
    {500, 600, 500, 6.967e-13, 63.822e0},
    {600, 700, 600, 1.454e-13, 71.835e0},
    {700, 800, 700, 3.614e-14, 88.667e0},
    {800, 900, 800, 1.170e-14, 124.64e0},
    {900, 1000, 900, 5.245e-15, 181.05e0},
    {1000, -1, 1000, 3.019e-15, 268.00e0},
};

const int table_sz = sizeof(AtmTable) / sizeof(AtmLayer);

double
dso::air_density_models::exponential::density(double sat_altitude_km) noexcept {
  assert(sat_altitude_km>0e0);

  // compute altitude, found by subtracting the Earth’s radius from the
  // satellite’s radius
  const int alt = std::trunc(sat_altitude_km);
  int layer_nr = table_sz - 1;
  for (int i = 0; i < table_sz - 1; i++) {
    if (alt >= AtmTable[i].limit_low && alt < AtmTable[i].limit_high) {
      layer_nr = i;
      break;
    }
  }
  return AtmTable[layer_nr].nominal_density *
         std::exp(-(sat_altitude_km -
                    static_cast<double>(AtmTable[layer_nr].base_altitude)) /
                  AtmTable[layer_nr].scale_height); // [kg/m3]
}
