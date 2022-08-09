#include "eop.hpp"
#include "iers2010/iers2010.hpp"

int dso::EopLookUpTable::interpolate(double fmjd_utc, double &xp, double &yp,
                                     double &dut1) const noexcept {

  // apply corrections; first use INTERP to interpolate x,y,ut1 values for
  // the given date (in [mas], [msec])
  if (iers2010::interp_pole(mjd, xpa, ypa, ut1a, sz, fmjd_utc, xp, yp, dut1))
    return 1;
  
  // account for variations in polar motion (Dx,Dy) ocean-tides; results in
  // [μas] and [μsec]
  double dxp, dyp, dut2;
  if (iers2010::ortho_eop(fmjd_utc, dxp, dyp, dut2))
    return 2;

  // transform to [mas] [msec]
  xp += dxp * 1e-3;
  yp += dyp * 1e-3;
  dut1 += dut2 * 1e-3;

  // account for libration effects; results in [μas]
  if (iers2010::pmsdnut2(fmjd_utc, dxp, dyp))
    return 3;

  // transform to [mas]
  xp += dxp * 1e-3;
  yp += dyp * 1e-3;

  return 0;
}
