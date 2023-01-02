#include "datetime/dtcalendar.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iau.hpp"
#include <cmath>
#include <cstdio>
#include <datetime/dtfund.hpp>

double pms_gmst(double mjd0, double mjd1) noexcept {
  constexpr const double TWOPI = 6.283185307179586476925287e0;
  constexpr const double RAD2SEC = 86400e0 / TWOPI;
  constexpr const double rmjd0 = 51544.5e0;
  const double t = (mjd0 - rmjd0) / 36525e0 + (mjd1) / 36525e0;
  const double gmst =
      std::fmod(67310.54841e0 + t * ((8640184.812866e0 + 3155760000e0) +
                                     t * (0.093104e0 + t * (-0.0000062))),
                86400e0);
  return std::fmod(gmst / RAD2SEC, TWOPI);
}

double lib_gmst(double mjd0, double mjd1) noexcept {
  dso::TwoPartDate mjd(mjd0, mjd1);
  dso::TwoPartDate jd = mjd.jd_sofa();
  return iers2010::sofa::gmst06(jd._big, jd._small, jd._big, jd._small);
}

double eq5_gmst(double mjd0, double mjd1) noexcept {
  dso::TwoPartDate mjd(mjd0, mjd1);
  //dso::TwoPartDate jd = mjd.jd_sofa();
  
  const double tu = (mjd0 - dso::j2000_mjd) + mjd1;
  double ipart;
  const double utfraction = std::modf(mjd1+.5e0, &ipart);

  const double era = /* Eq. 5.15, [rad] */
      dso::anp<double>(dso::D2PI *
               (utfraction + 0.7790572732640e0 + 0.00273781191135448e0 * tu));
  // printf("Note: m_era=%.12f s_era=%.12f diff=%.3e\n", era, iers2010::sofa::era00(jd._big, jd._small), era-iers2010::sofa::era00(jd._big, jd._small));
  const double t = mjd.jcenturies_sinceJ2000();
  const double gmst =
      std::fmod(/* Eq. 5.32 [rad] */
                era +
                    dso::sec2rad(0.014506e0 +
                                 t * (4612.156534e0 +
                                      t * (1.3915817e0 +
                                           t * (-0.00000044e0 +
                                                t * (-0.000029956e0 +
                                                     t * (-0.0000000368e0)))))),
                dso::D2PI);
  return gmst;
}

double grp_gmst(double mjd0, double mjd1) noexcept {
  double Tu0 = (mjd0 - 51544.5) / 36525.0;

  double GMST0 = (6.0 / 24 + 41.0 / (24 * 60) + 50.54841 / (24 * 60 * 60)) +
                 (8640184.812866 / (24 * 60 * 60)) * Tu0 +
                 (0.093104 / (24 * 60 * 60)) * Tu0 * Tu0 +
                 (-6.2e-6 / (24 * 60 * 60)) * Tu0 * Tu0 * Tu0;
  double r = 1.002737909350795 + 5.9006e-11 * Tu0 - 5.9e-15 * Tu0 * Tu0;
  return fmod(dso::D2PI * (GMST0 + r * mjd1), dso::D2PI);
}

double i82_gmst(double mjd0, double mjd1) noexcept {
  /* Coefficients of IAU 1982 GMST-UT1 model */
  constexpr const double A = 24110.54841 - dso::sec_per_day / 2.0;
  constexpr const double B = 8640184.812866;
  constexpr const double C = 0.093104;
  constexpr const double D = -6.2e-6;

  dso::TwoPartDate mjd(mjd0, mjd1);
  const double t = mjd.jcenturies_sinceJ2000();

  /* Fractional part of JD(UT1), in seconds. */
  const double f =
      dso::sec_per_day * (std::fmod(mjd0 + .5e0, 1.0) + std::fmod(mjd1, 1.0));
  /* GMST at this UT1. */
  return dso::anp<double>(dso::sec2rad((A + (B + (C + D * t) * t) * t) + f));
}

int main() {
  // suppose we have some MJD in UT time:
  dso::TwoPartDate d(59941e0, 0e0);
  while (d._small < 1e0) {
    const double g1 = pms_gmst(d._big, d._small);
    const double g2 = grp_gmst(d._big, d._small);
    const double g3 = lib_gmst(d._big, d._small);
    const double g4 = eq5_gmst(d._big, d._small);
    const double g5 = i82_gmst(d._big, d._small);

    printf("%.9f %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n", d._small,
           dso::rad2sec(g1 - g2), dso::rad2sec(g1 - g3), dso::rad2sec(g1 - g4),
           dso::rad2sec(g1 - g5), dso::rad2sec(g2 - g3), dso::rad2sec(g2 - g4),
           dso::rad2sec(g3 - g4));
    d._small += 1e-3;
  }

  return 0;
}
