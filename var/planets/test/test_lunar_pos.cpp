#include "datetime/datetime_write.hpp"
#include "datetime/dtcalendar.hpp"
#include "planetpos.hpp"
#include <cassert>
#include <cstdio>

int main() {
  dso::datetime<dso::seconds> t(dso::year(2006), dso::month(3),
                                dso::day_of_month(14));
  const int n_steps = 24;
  double pos[3];
  printf(" Moon position from low precision analytical theory\n");
  printf(" Date [TT]                  Position [km]\n");

  auto t1 = t;
  for (int i = 0; i <= n_steps; i++) {
    dso::moon_vector(t1, pos);
    printf("%s %15.3f %15.3f %15.3f\n", dso::strftime_ymd_hms(t1).c_str(),
           pos[0]/1e3, pos[1]/1e3, pos[2]/1e3);
    t1.add_seconds(dso::seconds(86400L / 2L));
  }

  printf(" Sun  position from low precision analytical theory\n");
  printf(" Date [TT]                  Position [km]\n");
  t1 = t;
  for (int i = 0; i <= n_steps; i++) {
    dso::sun_vector(t1, pos);
    printf("%s %15.3f %15.3f %15.3f\n", dso::strftime_ymd_hms(t1).c_str(),
           pos[0]/1e3, pos[1]/1e3, pos[2]/1e3);
    t1.add_seconds(dso::seconds(86400L / 2L));
  }
  return 0;
}
