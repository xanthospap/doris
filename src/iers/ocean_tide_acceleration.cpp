#include "egravity.hpp"
#include "iers2010/iau.hpp"
#include "tides.hpp"

/* References:
 * [2] Wunsch 2008, Comparison of two different ocean tide models especially
 * with respect to the GRACE satellite mission https://d-nb.info/97498454x/34
 */

/*
 * @param[in] t_tt  datetime in TT
 * @param[in] ut1_mjd Corresponding datetime in UT1 time scale, as MJD
 * @param[in] rsat ECEF vector of satellite, as (X,Y,Z) in [m]
 * @param[out] acc Acceleration at given point, as (X,Y,Z) in [m/sec^2]
 * @param[out] gradient Gradient of acceleration (i.e. da/dr)
 * @param[in] max_degree Max degree for computation of geopotential 
 *            coefficients. If <0 (default), then we are using the instance's 
 *            max degree (stored in the dCS instance, i.e. dCS.max_degree())
 * @param[in] max_order Max order for computation of geopotential 
 *            coefficients. If <0 (default), then we are using the instance's 
 *            max order (stored in the dCS instance, i.e. dCS.max_order())
 * @param[in] V A dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> of
 *              size >= degree + 2, to compute SH. If not given, the function
 *              will allocate it.
 * @param[in] W A dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> of
 *              size >= degree + 2, to compute SH. If not given, the function
 *              will allocate it.
 */
dso::iStatus dso::OceanTide::acceleration(
    const dso::TwoPartDate &mjdtt, const dso::TwoPartDate &mjdut1,
    const Eigen::Matrix<double, 3, 1> &rsat, Eigen::Matrix<double, 3, 1> &acc,
    Eigen::Matrix<double, 3, 3> &gradient, int max_degree, int max_order,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *V,
    dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise> *W) noexcept {

  if (max_degree < 0)
    max_degree = this->max_degree();
  if (max_order < 0)
    max_order = this->max_order();

  /* check input degree and order */
  if (max_degree > this->max_degree() || max_order > this->max_order()) {
    fprintf(stderr,
            "[ERROR] Invalid degree/order for computation of ocean tidal "
            "geopotential (traceback: %s)\n",
            __func__);
    return dso::iStatus(1);
  }

  /* allocate workspace for SH computations if needed */
  if (!V) {
    V = new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
        max_degree + 3, max_degree + 3);
  }
  if (!W) {
    W = new dso::Mat2D<dso::MatrixStorageType::LwTriangularColWise>(
        max_degree + 3, max_degree + 3);
  }

  /* check sizes of workspace matrices */
  if (V->rows() < max_degree + 2 || W->rows() < max_degree + 2) {
    fprintf(stderr,
            "[ERROR] Invalid workspace matrices for SH computation "
            "(traceback: %s)\n",
            __func__);
    return dso::iStatus(1);
  }

  /* compute geopotential corrections ΔC and ΔS */
  if (this->operator()(mjdtt, mjdut1, max_degree, max_order))
    return dso::iStatus(1);

  /* compute acceleration at satellite position (ITRF, cartesian) */
  dso::gravity_acceleration(dCS, rsat, max_degree, dCS.Re(), dCS.GM(), acc,
                            gradient, V, W);

  return dso::iStatus(0);
}

/* @brief Compute corrections to normalized C and S gravitational
 *        coefficients, using the model(s) described in IERS2010 standards.
 *        The corrections are computed and set in the instance's cs (member)
 *        variable
 *
 * @param[in] t_tt  datetime in TT
 * @param[in] ut1_mjd Corresponding datetime in UT1 time scale, as MJD
 * @param[in] max_degree Max degree for computation of geopotential 
 *            coefficients. If <0 (default), then we are using the instance's 
 *            max degree (stored in the dCS instance, i.e. dCS.max_degree())
 * @param[in] max_order Max order for computation of geopotential 
 *            coefficients. If <0 (default), then we are using the instance's 
 *            max order (stored in the dCS instance, i.e. dCS.max_order())
 */
dso::iStatus dso::OceanTide::operator()(const dso::TwoPartDate &mjdtt,
                                        const dso::TwoPartDate &mjdut1,
                                        int max_degree,
                                        int max_order) noexcept {
  /* check given degree and order */
  if (max_degree < 0)
    max_degree = dCS.max_degree();
  if (max_order < 0)
    max_order = dCS.max_order();

  /* check degree and order */
  if (max_degree > dCS.max_degree() || max_order > dCS.max_order() ||
      max_degree < max_order) {
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
  for (int n = 0; n <= max_degree; n++)
    dCS.S(n, 0) = 0e0;

  return dso::iStatus(0);
}
