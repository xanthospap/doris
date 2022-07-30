#ifndef __DSO_NRLMSISE00_CVERS_HPP__
#define __DSO_NRLMSISE00_CVERS_HPP__

#include "nrlmsise00_const.hpp"
#include <cmath>
#include <cstdint>
#include <cstring>
#ifdef DEBUG
#include <cassert>
#endif

namespace dso {

namespace nrlmsise00 {

struct switches {
  static constexpr const int dim = 25;
  double sw[dim] = {0e0};
  double swc[dim] = {0e0};

  switches() noexcept { set_on(); }

  void set_switch(int index, int val) noexcept {
#ifdef DEBUG
    assert(index >= 0 && index < dim);
#endif
    sw[index] = val % 2;
    swc[index] = (val == 1 || val == 2) ? 1e0 : 0e0;
  }

  void set_null(int index = -1) noexcept {
    if (index < 0) {
      // std::memset(isw, 0, sizeof(int8_t) * dim);
      std::memset(sw, 0, sizeof(double) * dim);
      std::memset(swc, 0, sizeof(double) * dim);
    } else {
#ifdef DEBUG
      assert(index < dim);
#endif
      // isw[index] = 0;
      sw[index] = 0e0;
      swc[index] = 0e0;
    }
  }

  void set_on(int index = -1) noexcept {
    if (index < 0) {
      // for (int i=0; i<dim; i++) isw[i] = 1;
      for (int i = 0; i < dim; i++)
        sw[i] = 1e0;
      for (int i = 0; i < dim; i++)
        swc[i] = 1e0;
    } else {
#ifdef DEBUG
      assert(index < dim);
#endif
      // isw[index] = 1;
      sw[index] = 1e0;
      swc[index] = 1e0;
    }
  }
}; // dso::nrlmsise00::switches

struct ApArray {
  double a[7];
}; // dso::nrlmsise00::ApArray

void spline(const double *__restrict__ x, const double *__restrict__ y, int n,
            double yp1, double ypn, double *__restrict__ y2,
            double *work /*size >= n */) noexcept;
double splini(const double *__restrict__ xa, const double *__restrict__ ya,
              double *__restrict__ y2a, int n, double x) noexcept;

double splint(const double *__restrict__ xa, const double *__restrict__ ya,
              double *__restrict__ y2a, int n, double x) noexcept;

struct InParams {
  nrlmsise00::switches sw;
  int year{0};  // year, currently ignored 
  int doy;      // day of year 
  double sec;   // seconds in day (UT) 
  double alt;   // altitude in kilometers 
  double glat;  // geodetic latitude 
  double glon;  // geodetic longitude 
  double lst;   // local apparent solar time (hours), see note below 
  double f107A; // 81 day average of F10.7 flux (centered on doy) 
  double f107;  // daily F10.7 flux for previous day 
  double ap;    // magnetic index(daily) 
  bool meters_{false}; // request result in m^3 and kg/m^3
  nrlmsise00::ApArray aparr; // ap index array

  bool meters() const noexcept { return meters_; }
  void meters_on() noexcept { meters_ = true; }
  double fdoy() const noexcept { return (double)doy + sec / 86400e0; }
  void set_switches_off() noexcept { sw.set_null(); }
  void set_switches_on() noexcept { sw.set_on(); }
  void set_switch(int index, double val) noexcept { sw.set_switch(index, val); }
  void set_ap_array(const double *ap_array) noexcept {
    std::memcpy(this->aparr.a, ap_array, sizeof(double) * 7);
  }
}; // InParams

struct OutParams {
  double d[9]; /* densities */
  double t[2]; /* temperatures */
};             // OutParams

double ccor2(double alt, double r, double h1, double zh, double h2) noexcept;
double ccor(double alt, double r, double h1, double zh) noexcept;
double dnet(double &dd, double dm, double zhm, double xmm, double xm) noexcept;
inline double g0(double a, const double *p) noexcept {
  return (a - 4e0 +
          (p[25] - 1e0) * (a - 4e0 +
                           (std::exp(-std::abs(p[24]) * (a - 4e0)) - 1e0) /
                               std::abs(p[24])));
}
inline double sumex(double ex) noexcept {
  return 1e0 + (1e0 - std::pow(ex, 19e0)) / (1e0 - ex) * std::pow(ex, 0.5e0);
}
inline double sg0(double ex, const double *p, const double *ap) noexcept {
  return (g0(ap[1], p) + (g0(ap[2], p) * ex + g0(ap[3], p) * ex * ex +
                          g0(ap[4], p) * std::pow(ex, 3e0) +
                          (g0(ap[5], p) * std::pow(ex, 4e0) +
                           g0(ap[6], p) * std::pow(ex, 12e0)) *
                              (1e0 - std::pow(ex, 8e0)) / (1e0 - ex))) /
         sumex(ex);
}
} // namespace nrlmsise00

struct Nrlmsise00 {
  static const double ptm[10];
  static const double pdm[8][10];
  static const double pavgm[10];
  
  // two dim arrays
  double pma[10][100];
  double pd[9][150];
  double ptl[4][100];
  double plg[4][9];
  const double pdl[2][25];

  // one dim arrays ...
  double pt[150];
  double ps[150];
  double tn1[5];
  double tn2[4];
  double tn3[5];
  double tgn1[2];
  double tgn2[2];
  double tgn3[2];
  double apt[4];
  const double zn2[4] = {72.5e0, 55e0, 45e0, 32.5e0};
  const double zn3[5] = {32.5e0, 20e0, 15e0, 10e0, 0e0};
  // const double sam[100]; -> just a function

  double gsurf;
  double re;
  double dd;
  double dm28, tlb, xg0, dm28m;
  double dfa;
  double ctloc, stloc;
  double c2tloc, s2tloc;
  double s3tloc, c3tloc;
  double apdf;

  Nrlmsise00() noexcept;
  
  constexpr double sam(int i) const noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 100);
#endif
    return (i < 50);
  }

  double densm(double alt, double d0, double xm, double &tz) const noexcept;
  double densu(double alt, double dlb, double t1, double t2, double xm,
               double xalpha, double &tz, double zlb, double s2,
               const double *zn1) noexcept;
  double zeta(double zz, double zl) const noexcept {
    return (zz - zl) * (re + zl) / (re + zz);
  }
  double scalh(double alt, double xm, double temp) const noexcept {
    const double g = gsurf / std::pow(1e0 + alt / re, 2e0);
    return dso::nrlmsise00::r100gas * temp / (g * xm);
  }
  double glatf(double lat, double &gv) const noexcept;

  /// @brief Neutral Atmosphere Empirical Model from the surface to lower
  /// exosphere.
  /// 
  /// @param[in] in Includes:
  ///   * doy day of year
  ///   * sec seconds in day
  ///   * alt altitude above surface (km)
  ///   * g_lat geodetic latitude
  ///   * g_long geodetic longitude
  ///   * lst local apparent solar time (hours)
  ///   * f107A 81 day average of F10.7 flux (centered on doy)
  ///   * f107 daily F10.7 flux for previous day
  ///   * ap magnetic index array
  ///   * d density array
  ///   * t temperature array
  /// 
  /// The output variables are:
  ///  - out.d[0]: HE number density(cm-3)
  ///  - out.d[1]: O  number density(cm-3)
  ///  - out.d[2]: N2 number density(cm-3)
  ///  - out.d[3]: O2 number density(cm-3)
  ///  - out.d[4]: AR number density(cm-3)
  ///  - out.d[5]: total mass density(gm/cm3) [includes d[8] in td7d]
  ///  - out.d[6]: H Number density(cm-3)
  ///  - out.d[7]: N Number density(cm-3)
  ///  - out.d[8]: Anomalous oxygen number density(cm-3)
  ///  - out.t[0]: exospheric temperature
  ///  - out.t[1]: temperature at alt
  /// 
  /// @note d[5] is the sum of the mass densities of the
  ///        species labeled by indices 0-4 and 6-7 in output variable d.
  ///        This includes He, O, N2, O2, Ar, H, and N but does NOT include
  ///        anomalous oxygen (species index 8).
  int gtd7(const nrlmsise00::InParams *in, int mass,
           nrlmsise00::OutParams *out) noexcept;
  
  /// @brief Neutral Atmosphere Empirical Model from the surface to lower 
  /// exosphere, including the anomalous oxygen contribution.
  /// 
  ///  This subroutine provides Effective Total Mass Density for output
  ///  d[5] which includes contributions from "anomalous oxygen" which can
  ///  affect satellite drag above 500 km.
  int gtd7d(const nrlmsise00::InParams *in, int mass,
            nrlmsise00::OutParams *out) noexcept;
  
  double glob7s(const nrlmsise00::InParams *in, double *p) noexcept;
  
  int ghp7(const nrlmsise00::InParams *in, nrlmsise00::OutParams *out,
           double press) noexcept;
  
  int gts7(const nrlmsise00::InParams *in, int mass,
           nrlmsise00::OutParams *out) noexcept;
  
  double globe7(const nrlmsise00::InParams *in, double *p) noexcept;
}; // Nrlmsise00

} // namespace dso
#endif
