#include "datetime/dtcalendar.hpp"
#include "planetpos.hpp"
#include <cstdio>
#include <cassert>

int main() {
  dso::datetime<dso::seconds> t1(dso::year(2006), dso::month(3),
                                 dso::day_of_month(14));
  const double mjd0 = t1.as_mjd();
  const double n_steps = 8;
  double pos[3];
  assert(std::abs(mjd0 - 53808e0)<1e-15);

      for (int i = 0; i < n_steps; i++) {
    double mjd_tt = mjd0 + i * 0.5e0;
    dso::moon_vector(mjd_tt, pos);
    printf("Moon at MJD=%15.5f : %15.3f %15.3f %15.3f\n", mjd_tt, pos[0],
           pos[1], pos[2]);
  }

  return 0;
}