#include "occultation.hpp"

double dso::conic_shadow_factor(const Eigen::Matrix<double, 3, 1> &rsat,
                                const Eigen::Matrix<double, 3, 1> &rsun,
                                double Rs, double Re) noexcept {
  /* apparent radius of the Earth */
  const double te = std::asin(Re / rsat.norm());
  /* apparent radius of the Sun */
  const double ts = std::asin(Rs / (rsun - rsat).norm());
  /* apparent separation between the Earth and the Sun */
  const double tes = std::acos(-(rsat.transpose() * (rsun - rsat))(0,0) /
                               (rsat.norm() * (rsun - rsat).norm()));
  /* no occultation */
  if (te + ts <= tes) {
    return 1e0;
  } else if (te - ts > tes) {
    /* Sun circular disk inside Earth circular disk */
    return 0e0;
  } else if (ts - te > tes) {
    /* Earth circular disk inside Sun circular disk */
    return 1e0 - te * te / (ts * ts);
  } else {
    /* CAF */
    const double caf =
        std::acos((ts * ts + tes * tes - te * te) / (2e0 * ts * tes));
    /* CBD */
    const double cbd =
        std::acos((te * te + tes * tes - ts * ts) / (2e0 * te * tes));
    /* areas */
    const double Safc = .5e0 * (caf * ts * ts);
    const double Saec = .5e0 * (ts * std::sin(caf)) * (ts * std::cos(caf));
    const double Sbdc = .5e0 * (cbd * te * te);
    const double Sbec = 0.5e0 * (te * std::sin(cbd)) * (te * std::cos(cbd));
    const double S = 2e0 * (Safc - Saec) + 2e0 * (Sbdc - Sbec);
    /* return factor */
    return 1e0 - S / (iers2010::DPI * ts * ts);
  }
}
