#include "nrlmsise00.hpp"

using namespace dso::nrlmsise00;

/// @brief Turbopause correction for msis models
double dso::nrlmsise00::dnet(double &dd, double dm, double zhm, double xmm,
                             double xm) noexcept {
  const double a = zhm / (xmm - xm);
  double dnet;

  if (dm <= 0e0 || dd <= 0e0) {
    if (std::abs(dd) < nearzero && std::abs(dm) < nearzero)
      dd = 1e0;
    if (std::abs(dm) < nearzero) {
      dnet = dd;
    } else if (std::abs(dd) < nearzero) {
      dnet = dm;
    }
  } else {
    const double ylog = a * std::log(dm / dd);
    if (ylog < -10e0) {
      dnet = dd;
    } else if (ylog > 10e0) {
      dnet = dm;
    } else {
      dnet = dd * std::pow(1e0 + std::exp(ylog), 1e0 / a);
    }
  }

  return dnet;
}
