#ifndef __DSO_NRLMSISE00_CVERS_HPP__
#define __DSO_NRLMSISE00_CVERS_HPP__

#include "nrlmsise00_const.hpp"
#include <cmath>
#include <cstdint>
#ifdef DEBUG
#include <cassert>
#endif

namespace dso {

namespace nrlmsise00 {

struct switches {
  static constexpr const int dim = 25;

  int8_t isw[dim] = {0 /*, [0 ... dim - 1] = 1*/}; // aka the sav
  double sw[dim];
  double swc[dim];

  void tselec(const double *sv) noexcept {
    for (int i = 0; i < dim; i++) {
      isw[i] = (int)sv[i];
      sw[i] = std::fmod(sv[i], 2e0);
      if (std::abs(sv[i] - 1e0) < nearzero ||
          std::abs(sv[i] - 2e0) < nearzero) {
        swc[i] = 1e0;
      } else {
        swc[i] = 0e0;
      }
    }
  }

  bool operator==(const switches &other) const noexcept {
    for (int i = 0; i < dim; i++)
      if (isw[i] != other.isw[i])
        return false;
    for (int i = 0; i < dim; i++)
      if (sw[i] != other.sw[i])
        return false;
    for (int i = 0; i < dim; i++)
      if (swc[i] != other.swc[i])
        return false;
    return true;
  }

  bool operator!=(const switches &other) const noexcept {
    return !(this->operator==(other));
  }

  int8_t operator()(int i) const noexcept {
#ifdef DEBUG
    assert(i < dim);
#endif
    return isw[i];
  }
};

struct ApArray {
  double a[7];
};

void spline(const double *__restrict__ x, const double *__restrict__ y, int n,
            double yp1, double ypn, double *__restrict__ y2,
            double *work /*size >= n */) noexcept;
double splini(const double *__restrict__ xa, const double *__restrict__ ya,
              double *__restrict__ y2a, int n, double x) noexcept;

double splint(const double *__restrict__ xa, const double *__restrict__ ya,
              double *__restrict__ y2a, int n, double x) noexcept;

struct InParams {
  int year{0};  /* year, currently ignored */
  int doy;      /* day of year */
  double sec;   /* seconds in day (UT) */
  double alt;   /* altitude in kilometers */
  double glat;  /* geodetic latitude */
  double glon;  /* geodetic longitude */
  double lst;   /* local apparent solar time (hours), see note below */
  double f107A; /* 81 day average of F10.7 flux (centered on doy) */
  double f107;  /* daily F10.7 flux for previous day */
  double ap;    /* magnetic index(daily) */
  bool meters_{false};
  nrlmsise00::ApArray ap_a; /* see above */
  nrlmsise00::switches sw;

  bool meters() const noexcept { return meters_; }
  double fdoy() const noexcept { return (double)doy + sec / 86400e0; }

  bool is_equal(const InParams &p, double limit=nearzero) const noexcept {
    return (year == p.year && doy == p.doy && std::abs(sec - p.sec) < limit &&
            std::abs(alt - p.alt) < limit && std::abs(glat - p.glat) < limit &&
            std::abs(glon - p.glon) < limit && std::abs(lst - p.lst) < limit &&
            std::abs(f107A - p.f107A) < limit &&
            std::abs(f107 - p.f107) < limit && std::abs(ap - p.ap) < limit &&
            sw == p.sw);
  }

}; // InParams

struct OutParams {
  double d[9]; /* densities */
  double t[2]; /* temperatures */
};             // OutParams

double ccor2(double alt, double r, double h1, double zh, double h2) noexcept;
double ccor(double alt, double r, double h1, double zh) noexcept;
double dnet(double &dd, double dm, double zhm, double xmm, double xm) noexcept;
double g0(double a, const double *p) noexcept {
  return (a - 4e0 +
          (p[25] - 1e0) * (a - 4e0 +
                           (std::exp(-std::abs(p[24]) * (a - 4e0)) - 1e0) /
                               std::abs(p[24])));
}
double sumex(double ex) noexcept {
  return 1e0 + (1e0 - std::pow(ex, 19e0)) / (1e0 - ex) * std::pow(ex, 0.5e0);
}
double sg0(double ex, const double *p, const double *ap) noexcept {
  return (g0(ap[1], p) + (g0(ap[2], p) * ex + g0(ap[3], p) * ex * ex +
                          g0(ap[4], p) * std::pow(ex, 3e0) +
                          (g0(ap[5], p) * std::pow(ex, 4e0) +
                           g0(ap[6], p) * std::pow(ex, 12e0)) *
                              (1e0 - std::pow(ex, 8e0)) / (1e0 - ex))) /
         sumex(ex);
}
} // namespace nrlmsise00

struct Nrlmsise00 {
  /* STATIC */
  /* POWER7 */
static const double pt[150];
static const double pd[9][150];
static const double ps[150];
static const double pdl[2][25];
static const double ptl[4][100];
static const double pma[10][100];
static const double sam[100];

/* LOWER7 */
static const double ptm[10];
static const double pdm[8][10];
static const double pavgm[10];

/* PARMB */
  double gsurf;
  double re;
  /* GTS3C */
  double dd;
  /* DMIX */
  double dm04, dm16, dm28, dm32, dm40, dm01, dm14;
  /* MESO7 */
  double tn1[5];
  double tn2[4];
  double tn3[5];
  double tgn1[2];
  double tgn2[2];
  double tgn3[2];
  /* LPOLY */
  double dfa;
  double plg[4][9];
  double ctloc, stloc;
  double c2tloc, s2tloc;
  double s3tloc, c3tloc;
  double apdf, apt[4];
  // VAR
  double zn1[5] = {120e0, 110e0, 100e0, 90e0, 72.e0};
  const double zn2[4] = {72.5e0, 55e0, 45e0, 32e0};
  const double zn3[5] = {32.5e0, 20e0, 15e0, 10e0, 0e0};
  const int mt[11] = {48, 0, 4, 16, 28, 32, 40, 1, 49, 14, 17};
  const double altl[8] = {200e0, 300e0, 160e0, 250e0,
                          240e0, 450e0, 320e0, 450e0};
  const double alpha[9] = {-0.3e0, 0e0, 0e0, 0e0, 0.1e0, 0e0, -0.3e0, 0e0, 0e0};

  double densm(double alt, double d0, double xm, double &tz) const noexcept;
  double densu(double alt, double dlb, double t1, double t2, double xm,
               double xalpha, double &tz, double zlb, double s2) noexcept;
  double zeta(double zz, double zl) const noexcept {
    return (zz - zl) * (re + zl) / (re + zz);
  }
  double scalh(double alt, double xm, double temp) const noexcept {
    const double g = gsurf / std::pow(1e0 + alt / re, 2e0);
    return dso::nrlmsise00::r100gas * temp / (g * xm);
  }
  double glatf(double lat, double &gv) const noexcept;
  int gtd7(const nrlmsise00::InParams *in, int mass,
           nrlmsise00::OutParams *out) noexcept;
  int gtd7d(const nrlmsise00::InParams *in, int mass,
            nrlmsise00::OutParams *out) noexcept;
  double glob7s(const nrlmsise00::InParams *in, double *p) noexcept;
  int ghp7(const nrlmsise00::InParams *in, nrlmsise00::OutParams *out,
           double press) noexcept;
  int gts7(const nrlmsise00::InParams *in, int mass,
           nrlmsise00::OutParams *out) noexcept;
}; // Nrlmsise00

} // namespace dso
#endif
