#include "egravity.hpp"
#include "iers2010/iau.hpp"
#include "tides.hpp"

/* References:
 * [2] Wunsch 2008, Comparison of two different ocean tide models especially with 
 *     respect to the GRACE satellite mission
 *     https://d-nb.info/97498454x/34
 */

dso::iStatus dso::OceanTide::acceleration(
    const dso::TwoPartDate &mjdtt, const Eigen::Matrix<double, 3, 1> &rsat,
    Eigen::Matrix<double, 3, 1> &acc, Eigen::Matrix<double, 3, 3> &gradient,
    int max_degree, int max_order) noexcept {
  if (max_degree < 0)
    max_degree = this->max_degree();
  if (max_order < 0)
    max_order = this->max_order();

  /* compute geopotential corrections ΔC and ΔS */
  if (this->operator()(mjdtt, max_degree, max_order))
    return dso::iStatus(1);

  /* compute acceleration at satellite position (ITRF, cartesian) */
  dso::gravity_acceleration(dCS, rsat, max_degree, dCS.Re(), dCS.GM(), acc,
                            gradient, &V, &W);

  return dso::iStatus(0);
}

dso::iStatus dso::OceanTide::operator()(const dso::TwoPartDate &mjdtt,
                                        int max_degree,
                                        int max_order) noexcept {
  /* check given degree and order */
  if (max_degree < 0)
    max_degree = dCS.max_degree();
  if (max_order < 0)
    max_order = dCS.max_order();
  if (max_degree > dCS.max_degree() || max_order > dCS.max_order()) {
    fprintf(stderr,
            "[ERROR] Invalid degree/order for ocean tides computation "
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
      iers2010::sofa::fal03(mjdtt),
      /* mean anomaly of sun, l' */
      iers2010::sofa::falp03(mjdtt),
      /* F = L - Ω */
      iers2010::sofa::faf03(mjdtt),
      /* Mean Elongation of the Moon from the Sun, D */
      iers2010::sofa::fad03(mjdtt),
      /* Mean Longitude of the Ascending Node of the Moon, Ω */
      iers2010::sofa::faom03(mjdtt)};

  /* fundamental arguments to Doodson multipliers */
  double dmult[6];
  dso::fundarg2doodson(fundarg, dso::gmst_utc(mjdtt.tt2utc()), dmult);

  /* clear stokes coefficients to recomput */
  dCS.clear();

  /* iterate through main waves f and apply Eq. 6.15 for each f */
  for (const auto &f : doodsonFreqs) {
    /* angle θ(t) for wave f (i.e. f is DoodsonOceanTideConstituent) */
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

  /* Warning !!
   * According to IERS 2010, Sec. 6.3.2
   * Note that, for zonal terms, FES2004 takes the approach to set the 
   * retrograde coefficients C_f,n0^(-) and S_f,n0^(-) to zero and to double 
   * the prograde coefficients C_f,n0^(+) and S_f,n0^(+). Therefore, after 
   * applying Equation (6.15), the ΔC_n0 have the expected value but the ΔSn0 
   * must be set to zero.
   */
  // for (int n=0; n<=max_degree; n++) dCS.S(n,0) = 0e0;

  return dso::iStatus(0);
}

/* TODO What GMST anle should we use here? HARDISP and Groops(?) seem to
/// be using a different angle than the one defines in IERS2010
int dso::OceanTide::operator()(const dso::TwoPartDate &mjdtt, int max_degree,
                               int max_order) noexcept {
  // compute Julian centuries since J2000.0 (TT)
  const double t = mjdtt.jcenturies_sinceJ2000();

  // compute fundamental arguments (for given TT)
  const double fundarg[] = {
      iers2010::sofa::fal03(t),  // mean anomaly of moon, l
      iers2010::sofa::falp03(t), // mean anomaly of sun, l'
      iers2010::sofa::faf03(t),  // L - Ω, F
      iers2010::sofa::fad03(t),  // Mean Elongation of the Moon from the Sun, D
      iers2010::sofa::faom03(
          t) // Mean Longitude of the Ascending Node of the Moon, Ω
  };

  // compute GMST using IAU 2006/2000A [rad]
  const double gmst = dso::gmst_utc(mjdtt.tt2utc());
  // iers2010::sofa::gmst06(jdut._big, jdut._small, jdtt._big, jdtt._small);

  // Doodson fundamental arguments (β_i = [τ,s,h,p,N',pl])
  double beta[6];
  fundarg2doodson(fundarg, gmst, beta);

  // clear harmonics coefficients
  dCS.clear();

  // start iterating through all constituents ...
  for (const auto &ddson : doodsonFreqs) {
    // find angle θ_f(t) for given Doodson number
    const double theta = ddson.doodson_number().phase(beta);
    // trigonometric numbers
    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    // max degree and order
    const int maxDegree = std::min(max_degree, ddson.max_degree());
    const int maxOrder = std::min(max_order, ddson.max_order());
    // perform computation, column-wise
    for (int m = 0; m <= maxOrder; m++) {
      for (int n = m; n <= maxDegree; n++) {
        dCS.C(n, m) += ddson.delCp(n, m) * ct + ddson.delSp(n, m) * st +
                       ddson.delCm(n, m) * ct + ddson.delSm(n, m) * st;
        if (m) {
          dCS.S(n, m) += -ddson.delCp(n, m) * st + ddson.delSp(n, m) * ct +
                         ddson.delCm(n, m) * st - ddson.delSm(n, m) * ct;
        }
      }
    }
  }

  return 0;
}
*/
