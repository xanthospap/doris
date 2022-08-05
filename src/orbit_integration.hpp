#include "egravity.hpp"
#include "eop.hpp"
#include <cassert>

namespace dso {

struct IntegrationParameters {
  ///< time in TAI
  double mjd_tai;
  ///< EOP parameters Look-up table
  const dso::EopLookUpTable &eopLUT;
  ///< gravity harmonics
  const dso::HarmonicCoeffs &harmonics;
  ///< space for Lagrange polynomials
  dso::Mat2D<dso::MatrixStorageType::Trapezoid> *Lagrange_V{nullptr};
  dso::Mat2D<dso::MatrixStorageType::Trapezoid> *Lagrange_W{nullptr};
  ///< degree and order of geopotential harmonics
  int degree, order;

  IntegrationParameters(int degree_, int order_,
                        const dso::EopLookUpTable &eoptable_,
                        const dso::HarmonicCoeffs &harmonics_) noexcept
      : eopLUT(eoptable_), harmonics(harmonics_),
        Lagrange_V{new dso::Mat2D<dso::MatrixStorageType::Trapezoid>(
            degree_ + 3, order_ + 3)},
        Lagrange_W{new dso::Mat2D<dso::MatrixStorageType::Trapezoid>(
            degree_ + 3, order_ + 3)},
        degree(degree_), order(order_) {
    assert(degree_ == harmonics_.degree());
  };

  ~IntegrationParameters() noexcept {
    if (Lagrange_V)
      delete Lagrange_V;
    if (Lagrange_W)
      delete Lagrange_W;
  }
}; // Integration Parameters
} // namespace dso
