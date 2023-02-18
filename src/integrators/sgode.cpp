#include "sgode.hpp"

dso::SGOde::SGOde(dso::ODEfun _f, int _neqn, double rerr, double aerr,
                  dso::IntegrationParameters *_params) noexcept
    : f(_f), neqn(_neqn), iflag(dso::SGOde::IFLAG::RESTART), relerr(rerr),
      abserr(aerr), params(_params) {
  Phi = Eigen::MatrixXd(neqn, 16);
  ArraysNeqn = Eigen::MatrixXd(neqn, 5); // wt, p, ypout, yp, yy
  Arrays13 = new double[13 * 7]; // psi, alpha, beta, v, w, sig, g (col-major)
}

dso::SGOde::~SGOde() noexcept {
  if (Arrays13)
    delete[] Arrays13;
}
