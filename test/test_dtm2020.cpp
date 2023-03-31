#include "astrodynamics.hpp"
#include "beacon_tbl.hpp"
#include "datetime/datetime_write.hpp"
#include "datetime/utcdates.hpp"
#include "doris_observation_equations.hpp"
#include "doris_rinex.hpp"
#include "doris_utils.hpp"
#include "filters/filters.hpp"
#include "geodesy/geoconst.hpp"
#include "geodesy/units.hpp"
#include "iers2010/cel2ter.hpp"
#include "iers2010/hardisp.hpp"
#include "iers2010/iers2010.hpp"
#include "iers2010/iersc.hpp"
#include "iers2010/tropo.hpp"
#include "integrators.hpp"
#include "orbit_integration.hpp"
#include "satellites/jason3.hpp"
#include "satellites/jason3_quaternions.hpp"
#include "sp3/sp3.hpp"
#include "sp3/sv_interpolate.hpp" /* debug mode */
#include "tides.hpp"
#include "var_utils.hpp"
#include <cassert>
#include <cstdio>
#include <cstring>
#include <datetime/dtcalendar.hpp>
#include <datetime/dtfund.hpp>
#include <iers2010/eop.hpp>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [DTM_2020_F107_Kp.dat]\n", argv[0]);
    return 1;
  }

  dso::Dtm2020 dtm(argv[1]);
  double f[] = {80e0,0e0, 180e0,0e0}; // fl, fh
  double fbar[] = {80e0,0e0, 180e0,0e0};
  double akp[] = {3e0,0e0,3e0,0e0};
  dtm.debug_set(0e0,0e0,300e0,3.1415e0,180e0,f,fbar,akp);
  dtm.dtm3();
  printf("Low  %.3f %.3f %.5e\n", dtm.temperature(), dtm.exospheric_temperature(), dtm.total_density_grcm3());
  dtm.debug_set(0e0,0e0,300e0,3.1415e0,180e0,f+2,fbar+2,akp);
  dtm.dtm3();
  printf("High %.3f %.3f %.5e\n", dtm.temperature(), dtm.exospheric_temperature(), dtm.total_density_grcm3());
  return 0;
}
