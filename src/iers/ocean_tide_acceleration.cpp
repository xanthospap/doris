#include "tides.hpp"
#include "egravity.hpp"
#include "iers2010/iau.hpp"

int dso::OceanTide::acceleration(const dso::TwoPartDate &mjdtt,
                                 /*const dso::TwoPartDate &mjdut1,*/
                                 const Eigen::Matrix<double, 3, 1> &rsat,
                                 Eigen::Matrix<double, 3, 1> &acc,
                                 int max_degree, int max_order) noexcept {
  if (max_degree<0) max_degree = this->max_degree();
  if (max_order<0) max_order = this->max_order();
  // compute geopotential corrections ΔC and ΔS
  this->operator()(mjdtt,/*mjdut1,*/max_degree,max_order);
  // dCS.scale(1e-11);

  //printf("DC coefficients:\n");
  //printf("DC(0,0) = %.6e\n", dCS.C(0,0));
  //printf("DC(1,0) = %.6e DC(1,1) = %.6e\n", dCS.C(1,0), dCS.C(1,1));
  //printf("DC(2,0) = %.6e DC(2,1) = %.6e DC(2,2) = %.6e\n", dCS.C(2,0), dCS.C(2,1), dCS.C(2,2));

  // compute acceleration at satellite position (ITRF, cartesian)
  Eigen::Matrix<double, 3, 3> partials;
  test::gravacc3(dCS, rsat, max_degree, dCS._Re, dCS._GM, acc, partials, &V, &W);

  return 0;
}

/// TODO What GMST anle should we use here? HARDISP and Groops(?) seem to 
/// be using a different angle than the one defines in IERS2010
int dso::OceanTide::operator()(const dso::TwoPartDate &mjdtt,
                               /*const dso::TwoPartDate &mjdut1,*/
                               int max_degree, int max_order) noexcept {
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
  // const auto jdtt = mjdtt.jd_sofa();
  // const auto jdut = mjdut1.jd_sofa();
  const double gmst =
      /* iers2010::sofa::gmst06(jdut._big, jdut._small, jdtt._big, jdtt._small); */
      dso::gmst_utc(mjdtt.tt2utc());

  // Doodson fundamental arguments (β_i = [τ,s,h,p,N',pl])
  double beta[6];
  fundarg2doodson(fundarg, gmst, beta);

  // clear harmonics coefficients
  dCS.clear();

  // start iterating through all constituents ...
  for (const auto &ddson : doodsonFreqs) {
    // find angle θ_f(t) for given Doodson number
    const double theta = ddson.doodson_number()./*theta_angle*/phase(beta);
    // trigonometric numbers
    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    // max degree and order
    const int maxDegree = std::min(max_degree, ddson.max_degree());
    const int maxOrder  = std::min(max_order, ddson.max_order());
    // perform computation, column-wise
    for (int m=0; m<=maxOrder; m++) {
      for (int n=m; n<=maxDegree; n++) {
        // printf("DC+(%d,%d) = %.5e\n", n,m,ddson.delCp(n, m));
        //if (n==0 && m==0) printf("DC(%d,%d) = %.6e + %.6e*%.3f + %.6e*%.3f + %.6e*%.3f + %.6e*%.3f\n", n,m,dCS.C(n, m),ddson.delCp(n, m), ct, ddson.delSp(n, m),st,ddson.delCm(n, m), ct, ddson.delSm(n, m), st);
        dCS.C(n, m) += ddson.delCp(n, m) * ct + ddson.delSp(n, m) * st +
                       ddson.delCm(n, m) * ct + ddson.delSm(n, m) * st;
        if (m) {
          dCS.S(n, m) += -ddson.delCp(n, m) * st + ddson.delSp(n, m) * ct +
                         ddson.delCm(n, m) * st - ddson.delSm(n, m) * ct;
          //if (n==2 && m==1) printf("DS(%d,%d) = %.6e\n", n,m,dCS.S(n, m));
        }
      }
    }
  }

  return 0;
}
