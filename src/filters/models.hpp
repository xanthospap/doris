#ifndef __DEFINE_MATH_MODEL_TYPES__
#define __DEFINE_MATH_MODEL_TYPES__

#include <cstring>
#include <datetime/dtcalendar.hpp>
#include <type_traits>

namespace dso {

template <typename XType> struct LinearModel {
  XType xref;
  double a{0e0}, b{0e0};

  int order() const noexcept { return 1; }
  double value_at(const XType &x) const noexcept { return a + b * (x - xref); }
}; // LinearModel

template <typename XType> struct PolynomialModel {
  int _order;
  XType xref;
  double *cf{nullptr};

  void resize(int new_order) noexcept {
    if (new_order != _order) {
      delete[] cf;
      cf = new double[new_order + 1];
      _order = new_order;
    }
  }

  inline double deltax(const XType &x1, const XType &x2) const noexcept {
    if constexpr (std::is_arithmetic_v<XType>)
      return x1 - x2;
    else
      return dso::date_diff<dso::DateTimeDifferenceType::FractionalSeconds>(x1,
                                                                            x2);
  }

  PolynomialModel(int porder) noexcept
      : _order(porder), cf(new double[porder + 1]){};

  ~PolynomialModel() noexcept {
    if (_order)
      delete[] cf;
  }

  PolynomialModel(const PolynomialModel &p) noexcept
      : _order(p._order), xref(p.xref), cf(new double[p._order + 1]) {
    std::memcpy(cf, p.cf, sizeof(double) * (_order + 1));
  }

  PolynomialModel &operator=(const PolynomialModel &p) noexcept {
    if (*this != p) {
      this->resize(p._order);
      std::memcpy(cf, p.cf, sizeof(double) * (_order + 1));
      xref = p.xref();
    }
    return *this;
  }

  PolynomialModel(PolynomialModel &&p) noexcept
      : _order(p._order), xref(p.xref), cf(p.cf){};

  PolynomialModel &operator=(PolynomialModel &&p) noexcept {
    _order = p._order;
    xref = p.xref;
    cf = p.cf;
    return *this;
  }

  int order() const noexcept { return _order; }

  double value_at(const XType &x) const noexcept {
    const double t = deltax(x, xref); // x - xref;
    double y = cf[0];
    double tpower = t;
    for (int i = 1; i < _order + 1; i++) {
      y += cf[i] * tpower;
      tpower *= t;
    }
    return y;
  }

}; // LinearModel
} // namespace dso
#endif
