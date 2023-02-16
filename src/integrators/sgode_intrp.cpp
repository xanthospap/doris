#include "sgode.hpp"
#ifdef DEBUG
#include <cassert>
#endif

///< According to Shampine & Gordon, intrp is an interpolation code for 
///< obtaining the solution at specified output points

/// x = member x
/// y = member yy
int dso::SGOde::intrp(double xout, Eigen::VectorXd &yout/*,
                      Eigen::Ref<Eigen::VectorXd> ypout*/) noexcept {
  double MemPool[3 * 13];
  double *__restrict__ g = MemPool;
  double *__restrict__ w = MemPool + 13;
  double *__restrict__ rho = MemPool + 26;

  const double hi = xout - x;
  const int ki = kold + 1;
  const int kip1 = ki + 1;

  g[0] = 1e0;
  rho[0] = 1e0;

  // initialize w(*) for computing g(*)
  for (int i = 0; i < ki; i++)
    w[i] = 1e0 / (double)(i + 1);

  // compute g(*)
  double term = 0e0;
  for (int j = 1; j < ki; j++) {
    const int jm1 = j - 1;
    const double psijm1 = psi(jm1);
    const double gamma = (hi + term) / psijm1;
    const double eta = hi / psijm1;
    const int limit = kip1 - j;
    for (int i = 0; i < limit - 1; i++) {
      w[i] = gamma * w[i] - eta * w[i + 1];
    }
    g[j] = w[0];
    rho[j] = gamma * rho[jm1];
    term = psijm1;
  }

  // interpolate
  ypout().setZero();
  yout.setZero();
  for (int j = 0; j < ki; j++) {
    const int i = kip1 - j - 2;
    const double t2 = g[i];
    const double t3 = rho[i];
#ifdef DEBUG
    assert(i < 16);
#endif
    yout += (t2 * Phi.col(i));
    ypout() += (t3 * Phi.col(i));
  }

  yout = yy() + hi * yout;

  return 0;
}
