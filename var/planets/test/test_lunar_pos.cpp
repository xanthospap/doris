#include "datetime/datetime_write.hpp"
#include "datetime/dtcalendar.hpp"
#include "planetpos.hpp"
#include <cassert>
#include <cstdio>

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s [DE/SPK KERNEL] [LEAPSEC/LSK KERNEL]\n",
            argv[0]);
    return 1;
  }
  dso::datetime<dso::seconds> t(dso::year(2006), dso::month(3),
                                dso::day_of_month(14));
  const int n_steps = 24;
  double pos[3], cpos[3];

  // to compute the planet's position via cspice, we need to load:
  // 1. the planetary ephemeris (SPK) kernel
  // 2. the leap-second (aka LSK) kernel
  dso::cspice::load_if_unloaded_spk(argv[1]);
  dso::cspice::load_if_unloaded_lsk(argv[2]);

  printf(" Moon position from low precision analytical theory\n");
  printf(" Date [TT]                  Position [km]\n");

  auto t1 = t;
  for (int i = 0; i <= n_steps; i++) {
    dso::moon_vector_approx(t1, pos);
    dso::moon_vector_cspice(t1, cpos);
    printf(
        "%s %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f\n",
        dso::strftime_ymd_hms(t1).c_str(), cpos[0], cpos[1], cpos[2], pos[0],
        pos[1], pos[2], cpos[0] - pos[0], cpos[1] - pos[1], cpos[2] - pos[2]);

    t1.add_seconds(dso::seconds(86400L / 2L));
  }

  printf(" Sun  position from low precision analytical theory\n");
  printf(" Date [TT]                  Position [km]\n");
  t1 = t;
  for (int i = 0; i <= n_steps; i++) {
    dso::sun_vector_approx2(t1, pos);
    dso::sun_vector_cspice(t1, cpos);
    printf(
        "%s %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f %15.3f\n",
        dso::strftime_ymd_hms(t1).c_str(), cpos[0], cpos[1], cpos[2], pos[0],
        pos[1], pos[2], cpos[0] - pos[0], cpos[1] - pos[1], cpos[2] - pos[2]);

    t1.add_seconds(dso::seconds(86400L / 2L));
  }
  return 0;
}
