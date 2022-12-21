#include "datetime/dtcalendar.hpp"
#include "datetime/utcdates.hpp"
#include "eop.hpp"
#include "iers2010/iers2010.hpp"
#include <algorithm>

namespace {
int lag_int(const dso::TwoPartDate &fmjd,
            const std::vector<dso::TwoPartDate> &t,
            const std::vector<double> &y, int order, int &index,
            double &val) noexcept {
  constexpr const dso::DateTimeDifferenceType FD =
      dso::DateTimeDifferenceType::FractionalDays;

  // find suitable index in input dates
  if (index < 0) {
    auto it = std::lower_bound(t.begin(), t.end(), fmjd);
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
      pval *= fmjd.diff<FD>(t[j]) / t[i].diff<FD>(t[j]);
    }
    for (int j = i + 1; j < index + window; j++) {
      pval *= fmjd.diff<FD>(t[j]) / t[i].diff<FD>(t[j]);
    }
    val += l * pval;
  }

  return 0;
}
} // unnamed namespace

int dso::EopLookUpTable::interpolate_lagrange(const dso::TwoPartDate &fmjd,
                                              dso::EopRecord &eopr,
                                              int order) const noexcept {
  int status = 0;
  int index = -1;

  // xp (pole)
  status += lag_int(fmjd, t, xp, order, index, eopr.xp);
  // yp (pole)
  status += lag_int(fmjd, t, yp, order, index, eopr.yp);
  // Dut1
  status += lag_int(fmjd, t, dut1, order, index, eopr.dut);
  // dX
  status += lag_int(fmjd, t, dX, order, index, eopr.dx);
  // dY
  status += lag_int(fmjd, t, dY, order, index, eopr.dy);
  // LOD
  status += lag_int(fmjd, t, lod, order, index, eopr.lod);

  // remember to assign date to the filled-in instance
  eopr.mjd = fmjd;

  return 0;
}

// do not regularize EOPs before calling this function
int dso::EopLookUpTable::interpolate(const dso::TwoPartDate &fmjd,
                                     dso::EopRecord &eopr,
                                     int order) const noexcept {
  
  // perform simple interpolation (lagrangian)
  dso::EopLookUpTable::interpolate_lagrange(fmjd, eopr, order);

//{
//  printf("[DEBUG@] utc=%.15f tt =%.15f\n", utc_fmjd.mjd(), tt_fmjd.mjd());
//  printf("[interp] xp = %.15f asec\n", eopr.xp);
//  printf("[interp] yp = %.15f asec\n", eopr.yp);
//  printf("[interp] dut= %.15f sec\n", eopr.dut);
//}

  // call ortho_eop
  double dxoc, dyoc, dut1oc;
  iers2010::ortho_eop(fmjd, dxoc, dyoc, dut1oc); // [μas] and [μsec]
//{
//    printf("\t[ortho_eop] xp =%.15f\n", dxoc);
//    printf("\t[ortho_eop] yp =%.15f\n", dyoc);
//    printf("\t[ortho_eop] dUT=%.15f\n", dut1oc);
//}

  // compute fundamental arguments (and gmst+π) needed for pmsdnut2 and
  // utlibr
  double fargs[6];
  iers2010::utils::eop_fundarg(fmjd, fargs);

  double dxlib, dylib;
  // iers2010::pmsdnut2(tt_fmjd, dxlib, dylib); // [μas]
  iers2010::utils::pmsdnut2(fmjd, fargs, dxlib, dylib);
//  {
//    printf("\t[pmsdnut2] xp =%.15f\n", dxlib);
//    printf("\t[pmsdnut2] yp =%.15f\n", dylib);
//  }

  double dut1lib, dlodlib;
  // iers2010::utlibr(tt_fmjd, dut1lib, dlodlib);
  iers2010::utils::utlibr(fmjd, fargs, dut1lib, dlodlib);
//{
//  printf("\t[utlibr] xp =%.15f\n", dut1lib);
//  printf("\t[utlibr] yp =%.15f\n", dlodlib);
//}

  // corrections in [asec]
  const double dx = (dxoc + dxlib) * 1e-6;
  const double dy = (dyoc + dylib) * 1e-6;
  const double dut = (dut1oc + dut1lib) * 1e-6;
  const double dlod = (dlodlib)*1e-6;

  // add corrections
  eopr.xp += dx;
  eopr.yp += dy;
  eopr.dut += dut;
  eopr.lod += dlod;

  return 0;
}
