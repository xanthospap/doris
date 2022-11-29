#include "datetime/utcdates.hpp"
#include "eop.hpp"
#include "iers2010/iers2010.hpp"
#include <datetime/dtfund.hpp>

int dso::EopLookUpTable::interpolate_lagrange(double tt_fmjd,
                                              dso::EopRecord &eopr,
                                              int order) const noexcept {
  int status = 0;
  int index = -1;

  // xp (pole)
  status +=
      iers2010::interp::lagint(mjd(), xp(), sz, tt_fmjd, eopr.xp, index, order);
  // yp (pole)
  status +=
      iers2010::interp::lagint(mjd(), yp(), sz, tt_fmjd, eopr.yp, index, order);
  // dUT1
  status += iers2010::interp::lagint(mjd(), dut(), sz, tt_fmjd, eopr.dut, index,
                                     order);
  // dX
  status +=
      iers2010::interp::lagint(mjd(), dx(), sz, tt_fmjd, eopr.dx, index, order);
  // dY
  status +=
      iers2010::interp::lagint(mjd(), dy(), sz, tt_fmjd, eopr.dy, index, order);
  // LOD
  status += iers2010::interp::lagint(mjd(), lod(), sz, tt_fmjd, eopr.lod, index,
                                     order);
  if (status) {
    fprintf(stderr,
            "[ERROR] Failed to interpolate EOP/ERP parameters for requested "
            "MJD=%.6f (traceback: %s)\n",
            tt_fmjd, __func__);
    return 1;
  }

  // remember to assign date to the filled-in instance
  eopr.mjd = tt_fmjd;

  return 0;
}

int dso::EopLookUpTable::interpolate(double tt_fmjd, dso::EopRecord &eopr,
                                     int order) const noexcept {

  // perform simple interpolation (lagrangian)
  dso::EopLookUpTable::interpolate_lagrange(tt_fmjd, eopr, order);

  // TT as Julian centuries since J2000 (try to keep accuracy)
  double imjd;
  const double tt_fday = std::modf(tt_fmjd, &imjd);
  const double jc2000 =
      (imjd + dso::mjd0_jd - dso::j2000_jd) / 36525e0 + tt_fday / 36525e0;

  // compute the effects of zonal Earth tides on the rotation of the Earth.
  // (assuming ERP values are regularized)
  double dut1, // [seconds]
      dlod,    // [seconds / day]
      domega;  // [radians / second]
  iers2010::rg_zont2(jc2000, dut1, dlod, domega);

  // add the effect to the ERP values
  eopr.dut += dut1;     // [seconds]
  eopr.lod += dlod;     // [seconds/day]
  eopr.omega += domega; // [?]

  // add effect of ocean tides (see iers2010::interp_pole and interp.f)
  double cx, cy, cdut, clod;
  iers2010::interp::pmut1_oceans(jc2000, cx, cy, cdut, clod);
  eopr.xp += cx;
  eopr.yp += cy;
  eopr.dut += cdut;
  eopr.lod += clod;

  // add lunisolar effect (libration)
  iers2010::interp::pm_gravi(jc2000, cx, cy);
  eopr.xp += cx;
  eopr.yp += cy;

  eopr.mjd = tt_fmjd; // TT MJD

  return 0;
}

void dso::EopLookUpTable::regularize() noexcept {
  dso::modified_julian_day tt_mjd;
  double imjd;
  for (int i = 0; i < sz; i++) {
    // TT as Julian centuries since J2000 (try to keep accuracy)
    const double tt_fday = std::modf(*mjd(i), &imjd);
    const double jc2000 =
        (imjd + dso::mjd0_jd - dso::j2000_jd) / 36525e0 + tt_fday / 36525e0;

    // compute the effects of zonal Earth tides on the rotation of the Earth.
    double dut1, // [seconds]
        dlod,    // [seconds / day]
        domega;  // [radians / second]
    iers2010::rg_zont2(jc2000, dut1, dlod, domega);

    // subtract the effect from the EOP values
    *(dut(i)) -= dut1;     // [seconds]
    *(lod(i)) -= dlod;     // [seconds/day]
    *(omega(i)) -= domega; // [?]
  }
  return;
}
