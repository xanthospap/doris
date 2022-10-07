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

int dso::EopLookUpTable::interpolate(double fmjd_utc, double &dx,
                                     double &dy) const noexcept {
  // quick return in case the passed in date is out of bounds ...
  if (fmjd_utc < *mjd || fmjd_utc >= mjd[sz - 1]) {
    fprintf(stderr,
            "[ERROR] Requested interpolation at %.5f but EOP values span [%.5f "
            "to %.5f] (traceback: %s)\n",
            fmjd_utc, mjd[0], mjd[sz - 1], __func__);
    return 1;
  }

  int status = 0;
  int index = -1;
  status += iers2010::interp::lagint(mjd, dxa, sz, fmjd_utc, dx, index);
  status += iers2010::interp::lagint(mjd, dya, sz, fmjd_utc, dy, index);
  if (status) {
    fprintf(stderr,
            "[ERROR] Failed EOP interpolation; requested mjd %.5f (traceback: "
            "%s)\n",
            fmjd_utc, __func__);
    return 9;
  }

  return 0;
}
