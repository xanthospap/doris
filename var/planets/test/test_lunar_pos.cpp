#include "datetime/datetime_write.hpp"
#include "datetime/dtcalendar.hpp"
#include "planetpos.hpp"
#include <cassert>
#include <cstdio>

double frac2sec(double mjd) noexcept {
  return (mjd - (int)mjd) * 86400;
}

int main() {
  dso::datetime<dso::seconds> t(dso::year(2006), dso::month(3),
                                 dso::day_of_month(14));
  const int n_steps = 8;
  double pos[3];

  for (int i = 0; i < n_steps; i++) {
    dso::moon_vector(t, pos);
    // printf("Moon at MJD=%15.5f : %15.3f %15.3f %15.3f\n", mjd_tt, pos[0],
    //       pos[1], pos[2]);
    printf("%s %15.5f %15.5f %15.5f\n", dso::strftime_ymd_hms(t).c_str(), pos[0], pos[1], pos[2]);
    t.add_seconds(dso::seconds(86400L/2L));
  }

  return 0;
}