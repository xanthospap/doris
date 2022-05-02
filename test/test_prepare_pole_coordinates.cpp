#include "iers2010/iers2010.hpp"
#include "iers_bulletin.hpp"
#include <datetime/dtcalendar.hpp>

// approximate number of data points in a bulletin B file (disregarding 
// preliminery values)
constexpr const int BBSZ = 40;

int main() {
  // construct a date ...
  dso::datetime<dso::milliseconds> t(dso::year{2022}, dso::month{1},
                                     dso::day_of_month{5},
                                     dso::milliseconds{3600L * 1000L});

  // get the corresponding bulletin b file from IERS (download)
  char bulletinb_fn[256];
  int error = dso::download_iers_bulletinb_for(t.mjd().as_underlying_type(),
                                               bulletinb_fn);
  if (error) {
    fprintf(stderr, "Failed to download Bulletin B file for date: %ld\n",
            t.mjd().as_underlying_type());
    return 1;
  }
  printf("Will try to compute the [X,Y] pole for mjd=%.5f\n", t.as_mjd());

  // parse the Bulletin B file (diregard preliminary values, only use final)
  dso::IersBulletinB_Section1Block bbblocks[35];
  dso::IersBulletinB bulb(bulletinb_fn);
  int bbblocks_size = bulb.parse_section1(bbblocks, false);
  if (bbblocks_size <= 0) {
    fprintf(stderr, "Failed parsing Bulletin B file %s\n", bulletinb_fn);
    return 2;
  }

  // apply corrections; first use INTERP to interpolate x,y,ut1 values for
  // the given date (in [mas], [msec])
  double xp, yp, ut1;
  double mjd[BBSZ], xpa[BBSZ], ypa[BBSZ], ut1a[BBSZ];
  for (int i = 0; i < bbblocks_size; i++)
    mjd[i] = static_cast<double>(bbblocks[i].mjd);
  for (int i = 0; i < bbblocks_size; i++)
    xpa[i] = bbblocks[i].x;
  for (int i = 0; i < bbblocks_size; i++)
    ypa[i] = bbblocks[i].y;
  for (int i = 0; i < bbblocks_size; i++)
    ut1a[i] = bbblocks[i].dut1;
  error = iers2010::interp_pole(mjd, xpa, ypa, ut1a, bbblocks_size, t.as_mjd(),
                                xp, yp, ut1);
  if (error) {
    fprintf(stderr,
            "Failed to interpolate EOP values for mjd=%.5f using Bulletin B "
            "file %s\n",
            t.as_mjd(), bulletinb_fn);
    return 3;
  }
  printf("%s %15.9f %15.9f %15.9f %15.3f %15.3f %15.3f [mas], [msec]\n",
         "interp    ", xp, yp, ut1, 0e0, 0e0, 0e0);

  // account for variations in polar motion (Dx,Dy) ocean-tides; results in
  // [μas]
  double dxp, dyp, dut1;
  if (iers2010::ortho_eop(t, dxp, dyp, dut1)) {
    fprintf(stderr, "Failed call to ortho_eop!\n");
    return 4;
  }
  xp += dxp * 1e-3;
  yp += dyp * 1e-3;
  ut1 += dut1 * 1e-3;
  printf("%s %15.9f %15.9f %15.9f %15.3f %15.3f %15.3f\n", "ortho_eop ", xp, yp,
         ut1, dxp * 1e-3, dyp * 1e-3, dut1 * 1e-3);

  // account for libration effects; results in [μas]
  if (iers2010::pmsdnut2(t, dxp, dyp)) {
    fprintf(stderr, "Failed call to pmsdnut2!\n");
    return 4;
  }
  xp += dxp * 1e-3;
  yp += dyp * 1e-3;
  printf("%s %15.9f %15.9f %15.9f %15.3f %15.3f %15.3f\n", "pmsdnut2  ", xp, yp,
         ut1, dxp * 1e-3, dyp * 1e-3, 0e0);

  // Final values stored in xp, yp in [mas]

  return 0;
}
