#include "eop.hpp"
#include <cstdio>
#include <datetime/dtfund.hpp>

int main(int argc, char *argv[]) {
  // need to provide a C04 file in cmds
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [EOP C04 file]\n", argv[0]);
    return 1;
  }

  // let's say we want EOP information for date:
  /*dso::datetime<dso::nanoseconds> t(dso::year(2022), dso::month(1),
                                    dso::day_of_month(1),
                                    dso::nanoseconds(60L * 1'000'000'000));*/
  dso::datetime<dso::nanoseconds> t(dso::year(2021), dso::month(1),
                                    dso::day_of_month(1),
                                    dso::nanoseconds(0));

  // extract information for +- 5 days (stop date is exclusive!)
  // dso::modified_julian_day mjd_start = t.mjd() - dso::modified_julian_day(5);
  // dso::modified_julian_day mjd_end = t.mjd() + dso::modified_julian_day(6);

  // ... or set some realy old years to plot results ...
  dso::datetime<dso::nanoseconds> tstart(dso::year(2021), dso::month(12),
                                    dso::day_of_month(1),
                                    dso::nanoseconds(0));
  dso::datetime<dso::nanoseconds> tend(dso::year(2022), dso::month(2),
                                    dso::day_of_month(1),
                                    dso::nanoseconds(0));
  // an EopLookUpTable to store results
  dso::EopLookUpTable eops;

  // parse the file and collect given dates ...
  if (dso::parse_iers_C04(argv[1], tstart.mjd(), tend.mjd(), eops)) {
    fprintf(stderr, "Failed parsing IERS/C04 file %s\n", argv[1]);
    return 1;
  }

  // regularize EOPs (DUT1 and DLOD)
  eops.regularize();

  // let's see what we have collected (just to verify)
  //for (int i = 0; i < eops.size(); i++) {
  //  printf("mjd %.1f xp %+10.8f yp %+10.8f dut1 %+10.8f lod %+10.8f dx %+10.8f "
  //         "dy %+10.8f\n",
  //         *eops.mjd(i), *eops.xp(i), *eops.yp(i), *eops.dut(i), *eops.lod(i),
  //         *eops.dx(i), *eops.dy(i));
  //}

  // interpolate for a range of 5 days, every 30 seconds
  dso::datetime<dso::nanoseconds> t_end(dso::year(2022), dso::month(1),
                                    dso::day_of_month(5),
                                    dso::nanoseconds(0));
  t = dso::datetime<dso::nanoseconds>(dso::year(2021), dso::month(12),
                                    dso::day_of_month(28),
                                    dso::nanoseconds(0));
  dso::EopRecord eopr;
  while (t < t_end) {
    if (eops.interpolate(t.as_mjd(), eopr)) {
      fprintf(stderr, "ERROR. Failed to interpolate at t=%.6f MJD\n", t.as_mjd());
      return 9;
    }
    printf("mjd %.9f xp %+10.8f yp %+10.8f dut1 %+10.8f lod %+10.8f dx %+10.8f "
           "dy %+10.8f\n",
           eopr.mjd, eopr.xp, eopr.yp, eopr.dut, eopr.lod, eopr.dx, eopr.dy);
    t.add_seconds(dso::seconds(30));
  }

  return 0;
}
