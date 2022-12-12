#include "datetime/utcdates.hpp"
#include "eop.hpp"
#include "iers2010/iers2010.hpp"
#include <datetime/dtfund.hpp>

int dso::EopLookUpTable::interpolate_lagrange(double tt_fmjd,
                                              dso::EopRecord &eopr,
                                              int order) const noexcept {
  int status = 0;
  int index = -1;
  double *tmpd = new double[order+1];

  // xp (pole)
  status +=
      iers2010::interp::lagint(mjd(), xp(), sz, tt_fmjd, eopr.xp, index, order);
  // yp (pole)
  status +=
      iers2010::interp::lagint(mjd(), yp(), sz, tt_fmjd, eopr.yp, index, order);

  int window = (order+1) / 2;
  int k = 0;
  for (int i=index-window; i<index+window; i++) {
    double it;
    double ft = std::modf(*mjd(i), &it);
    const double t = (it - dso::j2000_mjd) / dso::days_in_julian_cent +
                     ft / dso::days_in_julian_cent;
    double du,dl,dm;
    iers2010::rg_zont2(t, du, dl, dm);
    tmpd[k++] = *dut(i) - du;
  }
  // dUT1
  status += iers2010::interp::lagint(mjd(), tmpd-(index-window), sz, tt_fmjd, eopr.dut, index,
                                     order);
  {
  double it;
  double ft = std::modf(tt_fmjd, &it);
  const double t = (it - dso::j2000_mjd) / dso::days_in_julian_cent +
                   ft / dso::days_in_julian_cent;
  double du,dl,dm;
  iers2010::rg_zont2(t, du, dl, dm);
  eopr.dut += du;
  }
  
  //status += iers2010::interp::lagint(mjd(), dut(), sz, tt_fmjd, eopr.dut, index,
  //                                   order);
  // dX
  status +=
      iers2010::interp::lagint(mjd(), dx(), sz, tt_fmjd, eopr.dx, index, order);
  // dY
  status +=
      iers2010::interp::lagint(mjd(), dy(), sz, tt_fmjd, eopr.dy, index, order);
  // LOD
  status += iers2010::interp::lagint(mjd(), lod(), sz, tt_fmjd, eopr.lod, index,
                                     order);

  delete[] tmpd;
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

// this version is based on the interp.f routine but does not work well enough
int dso::EopLookUpTable::interpolate2(double tt_fmjd, dso::EopRecord &eopr,
                                     int order) const noexcept {

  // perform simple interpolation (lagrangian)
  dso::EopLookUpTable::interpolate_lagrange(tt_fmjd, eopr, order);

  // TT as Julian centuries since J2000 (try to keep accuracy)
  double it;
  double ft = std::modf(tt_fmjd, &it);
  const double t = (it - dso::j2000_mjd) / dso::days_in_julian_cent +
                   ft / dso::days_in_julian_cent;

  // compute the effects of zonal Earth tides on the rotation of the Earth.
  // (assuming ERP values are regularized)
  double dut1, // [seconds]
      dlod,    // [seconds / day]
      domega;  // [radians / second]
  iers2010::rg_zont2(t, dut1, dlod, domega);

  // add the effect to the ERP values
  eopr.dut += dut1;     // [seconds]
  eopr.lod += dlod;     // [seconds/day]
  eopr.omega += domega; // [?]

  // add effect of ocean tides (see iers2010::interp_pole and interp.f)
  double cx, cy, cdut, clod;
  iers2010::interp::pmut1_oceans(t, cx, cy, cdut, clod);
  eopr.xp += cx;
  eopr.yp += cy;
  eopr.dut += cdut;
  eopr.lod += clod;

  // add lunisolar effect (libration)
  iers2010::interp::pm_gravi(t, cx, cy);
  eopr.xp += cx;
  eopr.yp += cy;


  eopr.mjd = tt_fmjd; // TT MJD

  return 0;
}

// do not regularize EOPs before calling this function
int dso::EopLookUpTable::interpolate(double tt_fmjd, dso::EopRecord &eopr,
                                     int order) const noexcept {

  // perform simple interpolation (lagrangian)
  dso::EopLookUpTable::interpolate_lagrange(tt_fmjd, eopr, order);

  // call ortho_eop
  double dxoc,dyoc,dut1oc;
  iers2010::ortho_eop(tt_fmjd,dxoc,dyoc,dut1oc); // [μas] and [μsec]
  
  // compute fundamental arguments (and gmst+π) needed for pmsdnut2 and 
  // utlibr
  double fargs[6];
  iers2010::utils::eop_fundarg(tt_fmjd, fargs);

  double dxlib,dylib;
  //iers2010::pmsdnut2(tt_fmjd, dxlib, dylib); // [μas]
  iers2010::utils::pmsdnut2(tt_fmjd, fargs, dxlib, dylib);

  double dut1lib,dlodlib;
  // iers2010::utlibr(tt_fmjd, dut1lib, dlodlib);
  iers2010::utils::utlibr(tt_fmjd, fargs, dut1lib, dlodlib);

  // corrections in [asec]
  const double dx = (dxoc + dxlib) * 1e-6;
  const double dy = (dyoc + dylib) * 1e-6;
  const double dut = (dut1oc+dut1lib) * 1e-6;
  const double dlod = (dlodlib) * 1e-6;

  // add corrections
  eopr.xp += dx;
  eopr.yp += dy;
  eopr.dut += dut;
  eopr.lod += dlod;

  return 0;
}

void dso::EopLookUpTable::__regularize() noexcept {
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
