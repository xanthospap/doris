#include "eop.hpp"
#include "iers2010/iers2010.hpp"
#include "datetime/utcdates.hpp"
#include <datetime/dtfund.hpp>

int dso::EopLookUpTable::interpolate2(double fmjd_utc,
                                      dso::EopRecord &eopr) const noexcept {
  // interpolate regularized Dut1 values ...
  double corlod;
  int error = this->interpolate(fmjd_utc, eopr.xp, eopr.yp, eopr.ut1, corlod);

  double clod;
  int index = -1;
  if (iers2010::interp::lagint(mjd(), lod(), sz, fmjd_utc, clod, index)) {
    fprintf(stderr,
            "[ERROR] Failed EOP interpolation; requested mjd %.5f (traceback: "
            "%s)\n",
            fmjd_utc, __func__);
    return 9;
  }
  eopr.lod = clod + corlod;

  double futc_day, imjd, ftt_day;
  dso::modified_julian_day tt_mjd;
  // de-regularize dut1 and lod
  futc_day = std::modf(fmjd_utc, &imjd);
  ftt_day = dso::utc2tai(dso::modified_julian_day(imjd), futc_day, tt_mjd);
  ftt_day += 32184e-3 / 86400e0;
  // TT as Julian centuries since J2000
  const double jd =
      static_cast<double>(tt_mjd.as_underlying_type()) + dso::mjd0_jd;
  const double jc2000 = (jd - dso::j2000_jd) / 36525e0 + ftt_day / 36525e0;
  // compute the effects of zonal Earth tides on the rotation of the Earth.
  // Mind the units!!
  double dut1, // [seconds]
      dlod,    // [seconds / day]
      domega;  // [radians / second]
  iers2010::rg_zont2(jc2000, dut1, dlod, domega);
  // add the effect to the EOP values
  eopr.ut1 += dut1 * 1e-3; // [milliseconds]
  eopr.lod += dlod * 1e-3; // [milliseconds]

  error += interpolate(fmjd_utc, eopr.dx, eopr.dy);
  return error;
}

int dso::EopLookUpTable::interpolate(double fmjd_utc, double &xpole, double &ypole,
                                     double &dut1,
                                     double &corlod) const noexcept {

  // apply corrections; first use INTERP to interpolate x,y,ut1 values for
  // the given date (in [mas], [msec])
  if (iers2010::interp_pole(mjd(), xp(), yp(), dut(), sz, fmjd_utc, xpole,
                            ypole, dut1, corlod))
    return 1;
  
  // account for variations in polar motion (Dx,Dy) ocean-tides; results in
  // [μas] and [μsec]
  double dxp, dyp, dut2;
  if (iers2010::ortho_eop(fmjd_utc, dxp, dyp, dut2))
    return 2;

  // transform to [mas] [msec]
  xpole += dxp * 1e-3;
  ypole += dyp * 1e-3;
  dut1 += dut2 * 1e-3;

  // account for libration effects; results in [μas]
  if (iers2010::pmsdnut2(fmjd_utc, dxp, dyp))
    return 3;

  // transform to [mas]
  xpole += dxp * 1e-3;
  ypole += dyp * 1e-3;

  return 0;
}

int dso::EopLookUpTable::interpolate(double fmjd_utc, double &dxi,
                                     double &dyi) const noexcept {
  // quick return in case the passed in date is out of bounds ...
  if (fmjd_utc < *mjd() || fmjd_utc >= *mjd(sz - 1)) {
    fprintf(stderr,
            "[ERROR] Requested interpolation at %.5f but EOP values span [%.5f "
            "to %.5f] (traceback: %s)\n",
            fmjd_utc, *mjd(0), *mjd(sz - 1), __func__);
    return 1;
  }

  int status = 0;
  int index = -1;
  status += iers2010::interp::lagint(mjd(), dx(), sz, fmjd_utc, dxi, index);
  status += iers2010::interp::lagint(mjd(), dy(), sz, fmjd_utc, dyi, index);
  if (status) {
    fprintf(stderr,
            "[ERROR] Failed EOP interpolation; requested mjd %.5f (traceback: "
            "%s)\n",
            fmjd_utc, __func__);
    return 9;
  }

  return 0;
}

void dso::EopLookUpTable::regularize() noexcept {
  double futc_day, imjd, ftt_day;
  dso::modified_julian_day tt_mjd;
  
  for (int i=0; i<sz; i++) {
    // MJD is UTC at 0h 0min 0sec of day; transform to TT date
    futc_day = std::modf(*mjd(i), &imjd);
    ftt_day = dso::utc2tai(dso::modified_julian_day(imjd), futc_day, tt_mjd);
    ftt_day += 32184e-3 / 86400e0;

    // TT as Julian centuries since J2000
    const double jd =
        static_cast<double>(tt_mjd.as_underlying_type()) + dso::mjd0_jd;
    const double jc2000 = (jd - dso::j2000_jd) / 36525e0 + ftt_day / 36525e0;
    
    // compute the effects of zonal Earth tides on the rotation of the Earth.
    double dut1, // [seconds]
      dlod,      // [seconds / day]
      domega;    // [radians / second]
    iers2010::rg_zont2(jc2000, dut1, dlod, domega);

    // subtract the effect from the EOP values
    *(dut(i)) -= dut1; // [seconds]
    *(lod(i)) -= dlod; // [seconds/day]
    *(omega(i)) -= domega; // [?]
  }
  return;
}
