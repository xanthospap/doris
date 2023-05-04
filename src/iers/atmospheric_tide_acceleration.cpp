#include "egravity.hpp"
#include "iers2010/iau.hpp"
#include "tides.hpp"

/* References:
 * [2] Wunsch 2008, Comparison of two different ocean tide models especially
 * with respect to the GRACE satellite mission https://d-nb.info/97498454x/34
 */

dso::iStatus dso::AtmosphericTide::acceleration(
    const dso::TwoPartDate &mjdtt, const dso::TwoPartDate &mjdut1,
    const Eigen::Matrix<double, 3, 1> &rsat, Eigen::Matrix<double, 3, 1> &acc,
    Eigen::Matrix<double, 3, 3> &gradient, int max_degree,
    int max_order) noexcept {
  if (max_degree < 0)
    max_degree = this->max_degree();
  if (max_order < 0)
    max_order = this->max_order();

  /* compute geopotential corrections ΔC and ΔS */
  if (this->operator()(mjdtt, mjdut1, max_degree, max_order))
    return dso::iStatus(1);

  /* compute acceleration at satellite position (ITRF, cartesian) */
  dso::gravity_acceleration(dCS, rsat, max_degree, dCS.Re(), dCS.GM(), acc,
                            gradient, &V, &W);

  return dso::iStatus(0);
}

dso::iStatus dso::AtmosphericTide::operator()(const dso::TwoPartDate &mjdtt,
                                        const dso::TwoPartDate &mjdut1,
                                        int max_degree,
                                        int max_order) noexcept {
  /* check given degree and order */
  if (max_degree < 0)
    max_degree = dCS.max_degree();
  if (max_order < 0)
    max_order = dCS.max_order();
  if (max_degree > dCS.max_degree() || max_order > dCS.max_order()) {
    fprintf(stderr,
            "[ERROR] Invalid degree/order for atmospheric tides computation "
            "(traceback: %s)\n",
            __func__);
    return dso::iStatus(1);
  }

  /* we will need the the argument of the tide constituent θ(t) for every
   * wave f. Hence, compute the Doodson arguemnts for this epoch. Start with
   * fundamental arguemnts, fundarg
   */
  const double fundarg[] = {
      /* mean anomaly of moon, l */
      iers2010::fal03(mjdtt),
      /* mean anomaly of sun, l' */
      iers2010::falp03(mjdtt),
      /* F = L - Ω */
      iers2010::faf03(mjdtt),
      /* Mean Elongation of the Moon from the Sun, D */
      iers2010::fad03(mjdtt),
      /* Mean Longitude of the Ascending Node of the Moon, Ω */
      iers2010::faom03(mjdtt)};

  /* fundamental arguments to Doodson multipliers */
  double dmult[6];
  dso::fundarg2doodson(fundarg, iers2010::sofa::gmst06(mjdut1, mjdtt), dmult);

  /* clear stokes coefficients to recompute */
  dCS.clear();

  /* iterate through main waves f and apply Eq. 6.15 for each f */
  for (const auto &f : doodsonFreqs) {
    /* angle θ(t) for wave f (i.e. f is DoodsonAtmosphericTideConstituent) */
    const double theta = f.doodson_number().phase(dmult);
    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    /* iterate through n,m for current wave */
    for (int m = max_order; m >= 0; m--) {
      for (int n = max_degree; n >= m; n--) {
        /* ΔCnm,f = [C^(+)nm + C^(-)nm] * cosθ + [S^(+)nm + S^(-)nm] * sinθ
         * ΔSnm,f = [S^(+)nm - S^(-)nm] * cosθ - [C^(+)nm - C^(-)nm] * sinθ
         * see [2], Appendix A: Summation of the ocean tide height harmonics
         */
        dCS.C(n, m) += (ct * (f.delCp(n, m) + f.delCm(n, m)) +
                        st * (f.delSp(n, m) + f.delSm(n, m)));
        dCS.S(n, m) += (ct * (f.delSp(n, m) - f.delSm(n, m)) -
                        st * (f.delCp(n, m) - f.delCm(n, m)));
      } /* end looping degrees */
    }   /* end looping order */
  }     /* end current wave f */

  return dso::iStatus(0);
}
