#include "datetime/utcdates.hpp"
#include "eop.hpp"
#include "iers2010/iers2010.hpp"
#include "datetime/dtcalendar.hpp"
#include <algorithm>

namespace {
int lag_int(const dso::TwoPartDate &tt_fmjd,
            const std::vector<dso::TwoPartDate> &t,
            const std::vector<double> &y, int order, int &index,
            double &val) noexcept {
  constexpr const dso::DateTimeDifferenceType FD =
      dso::DateTimeDifferenceType::FractionalDays;

  // find suitable index in input dates
  if (index < 0) {
    auto it = std::lower_bound(t.begin(), t.end(), tt_fmjd);
    if (it == t.end())
      return 1;
    index = std::distance(t.begin(), it);
  }
  // make sure index is ok
  int window = (order + 1) / 2;
  if (index < window || index > (int)t.size() - window)
    return 1;
  // perform lagrangian interpolation, spanning indexes
  // [index-window, index+window)
  val = 0e0;
  for (int i = index - window; i < index + window; i++) {
    const double l = y[i];
    double pval = 1e0;
    for (int j = index - window; j < i; j++) {
      // pval *= (double)(tt_fmjd - t[j]) / (double)(t[i] - t[j]);
      pval *= tt_fmjd.diff<FD>(t[j]) / t[i].diff<FD>(t[j]);
    }
    for (int j = i + 1; j < index + window; j++) {
      // pval *= (double)(tt_fmjd - t[j]) / (double)(t[i] - t[j]);
      pval *= tt_fmjd.diff<FD>(t[j]) / t[i].diff<FD>(t[j]);
    }
    val += l * pval;
  }

  return 0;
}
}// unnamed namespace

int dso::EopLookUpTable::interpolate_lagrange(const dso::TwoPartDate &tt_fmjd,
                                              dso::EopRecord &eopr,
                                              int order) const noexcept {
  int status = 0;
  int index = -1;

  // xp (pole)
  status += lag_int(tt_fmjd, t, xp, order, index, eopr.xp);
  // yp (pole)
  status += lag_int(tt_fmjd, t, yp, order, index, eopr.yp);
  // Dut1
  status += lag_int(tt_fmjd, t, dut1, order, index, eopr.dut);
  // dX
  status += lag_int(tt_fmjd, t, dX, order, index, eopr.dx);
  // dY
  status += lag_int(tt_fmjd, t, dY, order, index, eopr.dy);
  // LOD
  status += lag_int(tt_fmjd, t, lod, order, index, eopr.lod);

  // remember to assign date to the filled-in instance
  eopr.mjd = tt_fmjd;

  return 0;
}

/* this version is based on the interp.f routine but does not work well enough
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
}*/

// do not regularize EOPs before calling this function
int dso::EopLookUpTable::interpolate(const dso::TwoPartDate &tt_fmjd,
                                     dso::EopRecord &eopr,
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

/*void dso::EopLookUpTable::__regularize() noexcept {
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
}*/
