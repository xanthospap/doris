#include "astrodynamics.hpp"
#include "datetime/datetime_write.hpp"
#include "egravity.hpp"
#include "eigen3/Eigen/Eigen"
#include "geodesy/units.hpp"
#include "icgemio.hpp"
#include "iers2010/iau.hpp"
#include "iers2010/iers2010.hpp"
#include "iers_bulletin.hpp"
#include "integrators.hpp"
#include "planetpos.hpp"
#include "sp3/sp3.hpp"
#include <cassert>
#include <cstdio>
#include <datetime/dtfund.hpp>

// usually using these datetimes ...
using Datetime = dso::datetime<dso::nanoseconds>;
using nsec_t = dso::nanoseconds::underlying_type;

int main(int argc, char *argv[]) {
     // handle gravity field
  if (argc != 3) {
    fprintf(stderr,
            "Usage: %s [SPK] [LSK]\n", argv[0]);
    return 1;
  }

  // to compute the planet's position via cspice, we need to load:
  // 1. the planetary ephemeris (SPK) kernel
  // 2. the leap-second (aka LSK) kernel
  dso::cspice::load_if_unloaded_spk(argv[1]);
  dso::cspice::load_if_unloaded_lsk(argv[2]);

  const Datetime mjd0 = Datetime(dso::year(2006), dso::month(3), dso::day_of_month(14),
                        dso::nanoseconds(0));


    printf("Date [TT]               Position [km]\n");
    double sat[3] = {-5174099.7245e0,   -3677755.6744e0,   +3661413.7036};

    const int    N_Step =   8;
    // const double Step   = 0.5; // [d]
    const nsec_t Step = 86400L / 2L; // half day in sec
    double r[3];
    Datetime mjd;
    char buf[64];
    for (int i=0; i<=N_Step; i++) {
        mjd = mjd0;
        mjd.add_seconds(dso::nanoseconds(static_cast<nsec_t>(i*Step * 1e9)));
        dso::cspice::j2planet_pos_from(dso::cspice::jd2et(mjd.as_jd()), 301, 399, r);
        const auto acc = dso::point_mass_accel(dso::Vector3::ptr2vec(sat), dso::Vector3::ptr2vec(r) * 1e3, 3.986004418e14);
        dso::strftime_ymd_hmfs(mjd, buf);
        printf("%s %+15.3f %+15.3f %+15.3f %+.6f %+.6f %+.6f\n", buf, r[0], r[1], r[2], acc(0), acc(1), acc(2));
    }

    return 0;
}