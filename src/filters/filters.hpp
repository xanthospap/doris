#ifndef __FILTER_IMPLEMENTATION_MODELS_HPP__
#define __FILTER_IMPLEMENTATION_MODELS_HPP__

#include "datetime/dtcalendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "filters/ekf.hpp"
#ifdef DEBUG
#include <cassert>
#include <cstdio>
#endif

namespace dso {

enum class BeaconClockModel : char { None, Linear };

template <int Np, BeaconClockModel bcm> struct EkfFilterImpl {
  ExtendedKalmanFilter _ekf;
  dso::TwoPartDate _ref_epoch;
  int _num_stations{0};

  static std::size_t num_params(int num_stations) noexcept {
    printf(
        ">> asked for num. params, answer is %d (computed with %d stations)\n",
        6 + Np + num_stations +
            num_stations * ((bcm == BeaconClockModel::None) ? (1) : (2)),
        num_stations);
    return 6              // satellite state vector
           + Np           // satellite force model parameters
           + num_stations // zenith tropospheric wet delay, per beacon
           // number of relative clock offset parameters, per beacon
           + num_stations * ((bcm == BeaconClockModel::None) ? (1) : (2));
  }

  int num_params() const noexcept {
    return EkfFilterImpl<Np, bcm>::num_params(_num_stations);
  }

  int num_stations() const noexcept { return _num_stations; }

  Eigen::MatrixXd estimates() const noexcept { return _ekf.x; }

  Eigen::MatrixXd &P() noexcept { return _ekf.P; }

  EkfFilterImpl(int num_stations,
                const dso::TwoPartDate &ref_epoch) noexcept
      : _ekf(EkfFilterImpl<Np, bcm>::num_params(num_stations)),
        _ref_epoch(ref_epoch), _num_stations(num_stations) {}

  void time_update(const dso::TwoPartDate &tk,
                   const Eigen::VectorXd &xk,
                   const Eigen::MatrixXd &phi) noexcept {
    _ekf.time_update(tk, xk, phi);
  }

  void observation_update(double z, double g, double sigma,
                          const Eigen::VectorXd &H) noexcept {
    _ekf.observation_update(z, g, sigma, H);
  }

}; // ekfFilterImpl

template <int Np, BeaconClockModel bcm> struct EkfFilter {};

template <int Np> struct EkfFilter<Np, BeaconClockModel::Linear> {
  EkfFilterImpl<Np, BeaconClockModel::Linear> _ekf;

  EkfFilter(int num_stations,
            const dso::TwoPartDate &ref_epoch) noexcept
      : _ekf(num_stations, ref_epoch) {}

  std::size_t num_params() const noexcept { return _ekf.num_params(); }

  Eigen::MatrixXd estimates() const noexcept { return _ekf.estimates(); }
  Eigen::MatrixXd &P() noexcept { return _ekf.P(); }

  double deltat(const dso::TwoPartDate &t) const noexcept {
    //const dso::nanoseconds dtsec = t.delta_sec(_ekf._ref_epoch);
    // return dtsec.to_fractional_seconds();
   return t.diff<dso::DateTimeDifferenceType::FractionalSeconds>(_ekf._ref_epoch);
  }

  constexpr std::size_t tropo_index(int beacon_nr) const noexcept {
    return 6 + Np + beacon_nr * 3;
  }

  constexpr std::size_t rfoff_index(int beacon_nr) const noexcept {
    return 6 + Np + beacon_nr * 3 + 1; // + dt * next-index !!
  }
  constexpr std::size_t drag_coef_index() const noexcept {
    static_assert(Np >= 1);
    return 6;
  }

  double drag_coef() const noexcept { return _ekf._ekf.x(drag_coef_index()); }
  double &drag_coef() noexcept { return _ekf._ekf.x(drag_coef_index()); }
  double tropo_estimate(int beacon_nr) const noexcept {
#ifdef DEBUG
    assert(tropo_index(beacon_nr) >= 0 &&
           tropo_index(beacon_nr) < _ekf._ekf.x.rows());
#endif
    return _ekf._ekf.x(tropo_index(beacon_nr));
  }
  double &tropo_estimate(int beacon_nr) noexcept {
    return _ekf._ekf.x(tropo_index(beacon_nr));
  }

  double
  rfoff_estimate(int beacon_nr,
                 const dso::TwoPartDate &t) const noexcept {
    const int idx = rfoff_index(beacon_nr);
#ifdef DEBUG
    assert(idx >= 0 && idx < _ekf._ekf.x.rows());
#endif
    //const dso::nanoseconds dtsec = t.delta_sec(_ekf._ref_epoch);
    //const double dt = dtsec.to_fractional_seconds();
    const double dt = t.diff<dso::DateTimeDifferenceType::FractionalSeconds>(_ekf._ref_epoch);
    return _ekf._ekf.x(idx) + dt * _ekf._ekf.x(idx + 1);
  }
  // double &rfoff_estimate(int beacon_nr) noexcept {
  //     return _ekf._ekf.x(rfoff_index(beacon_nr));
  // }

  void time_update(const dso::TwoPartDate &tk,
                   const Eigen::VectorXd &xk,
                   const Eigen::MatrixXd &phi) noexcept {
    _ekf.time_update(tk, xk, phi);
  }

  void observation_update(double z, double g, double sigma,
                          const Eigen::VectorXd &H) noexcept {
    _ekf.observation_update(z, g, sigma, H);
  }

  void set_tropo_apriori(double val, int beacon_nr = -1) noexcept {
#ifdef DEBUG
    assert(beacon_nr >= -1 && beacon_nr < _ekf._num_stations);
#endif
    if (beacon_nr > -1) {
      _ekf._ekf.x(tropo_index(beacon_nr)) = val;
    } else {
      const int start_idx = tropo_index(0);
      constexpr int step = 2;
      int idx = start_idx;
      for (int i = 0; i < _ekf.num_stations(); i++) {
        _ekf._ekf.x(idx) = val;
        idx += step;
      }
    }
  }

  void set_tropo_apriori_sigma(double val, int beacon_nr = -1) noexcept {
#ifdef DEBUG
    assert(beacon_nr >= -1 && beacon_nr < _ekf._num_stations);
#endif
    if (beacon_nr > -1) {
      _ekf._ekf.P(tropo_index(beacon_nr), tropo_index(beacon_nr)) = val;
    } else {
      const int start_idx = tropo_index(0);
      constexpr int step = 2;
      int idx = start_idx;
      for (int i = 0; i < _ekf.num_stations(); i++) {
        _ekf._ekf.P(idx, idx) = val;
        idx += step;
      }
    }
  }

  void set_rfoff_apriori(double a, double b, int beacon_nr = -1) noexcept {
#ifdef DEBUG
    assert(beacon_nr >= -1 && beacon_nr < _ekf._num_stations);
#endif
    if (beacon_nr > -1) {
      _ekf._ekf.x(rfoff_index(beacon_nr)) = a;
      _ekf._ekf.x(rfoff_index(beacon_nr) + 1) = b;

    } else {
      const int start_idx = rfoff_index(0);
      constexpr int step = 3;
      int idx = start_idx;
      for (int i = 0; i < _ekf.num_stations(); i++) {
        _ekf._ekf.x(idx) = a;
        _ekf._ekf.x(idx + 1) = b;
        idx += step;
      }
    }
  }

  void set_rfoff_apriori_sigma(double sa, double sb,
                               int beacon_nr = -1) noexcept {
#ifdef DEBUG
    assert(beacon_nr >= -1 && beacon_nr < _ekf._num_stations);
#endif
    if (beacon_nr > -1) {
      _ekf._ekf.P(rfoff_index(beacon_nr), rfoff_index(beacon_nr)) = sa;
      _ekf._ekf.P(rfoff_index(beacon_nr) + 1, rfoff_index(beacon_nr) + 1) = sb;
    } else {
      const int start_idx = rfoff_index(0);
      constexpr int step = 3;
      int idx = start_idx;
      for (int i = 0; i < _ekf.num_stations(); i++) {
        _ekf._ekf.P(idx, idx) = sa;
        _ekf._ekf.P(idx + 1, idx + 1) = sb;
        idx += step;
      }
    }
  }

  void set_drag_coef_apriori(double val) noexcept {
    _ekf._ekf.x(drag_coef_index()) = val;
  }
  void set_drag_coef_apriori_sigma(double val) noexcept {
    _ekf._ekf.P(drag_coef_index(), drag_coef_index()) = val;
  }
};

template <int Np> struct EkfFilter<Np, BeaconClockModel::None> {
  EkfFilterImpl<Np, BeaconClockModel::None> _ekf;

  EkfFilter(int num_stations,
            const dso::TwoPartDate &ref_epoch) noexcept
      : _ekf(num_stations, ref_epoch) {}

  std::size_t num_params() const noexcept { return _ekf.num_params(); }

  Eigen::MatrixXd estimates() const noexcept { return _ekf.estimates(); }
  Eigen::MatrixXd &P() noexcept { return _ekf.P(); }

  constexpr std::size_t tropo_index(int beacon_nr) const noexcept {
    return 6 + Np + beacon_nr * 2;
  }

  constexpr std::size_t rfoff_index(int beacon_nr) const noexcept {
    return 6 + Np + beacon_nr * 2 + 1;
  }
  constexpr std::size_t drag_coef_index() const noexcept {
    static_assert(Np >= 1);
    return 6;
  }

  double drag_coef() const noexcept { return _ekf._ekf.x(drag_coef_index()); }
  double &drag_coef() noexcept { return _ekf._ekf.x(drag_coef_index()); }
  double tropo_estimate(int beacon_nr) const noexcept {
    return _ekf._ekf.x(tropo_index(beacon_nr));
  }
  double &tropo_estimate(int beacon_nr) noexcept {
    return _ekf._ekf.x(tropo_index(beacon_nr));
  }
  double rfoff_estimate(int beacon_nr) const noexcept {
    return _ekf._ekf.x(rfoff_index(beacon_nr));
  }
  double &rfoff_estimate(int beacon_nr) noexcept {
    return _ekf._ekf.x(rfoff_index(beacon_nr));
  }

  void time_update(const dso::TwoPartDate &tk,
                   const Eigen::VectorXd &xk,
                   const Eigen::MatrixXd &phi) noexcept {
    _ekf.time_update(tk, xk, phi);
  }

  void observation_update(double z, double g, double sigma,
                          const Eigen::VectorXd &H) noexcept {
    _ekf.observation_update(z, g, sigma, H);
  }

  void set_tropo_apriori(double val, int beacon_nr = -1) noexcept {
#ifdef DEBUG
    assert(beacon_nr >= -1 && beacon_nr < _ekf._num_stations);
#endif
    if (beacon_nr > -1) {
      _ekf._ekf.x(tropo_index(beacon_nr)) = val;
    } else {
      const int start_idx = tropo_index(0);
      constexpr int step = 2;
      int idx = start_idx;
      for (int i = 0; i < _ekf.num_stations(); i++) {
        _ekf._ekf.x(idx) = val;
        idx += step;
      }
    }
  }

  void set_tropo_apriori_sigma(double val, int beacon_nr = -1) noexcept {
#ifdef DEBUG
    assert(beacon_nr >= -1 && beacon_nr < _ekf._num_stations);
#endif
    if (beacon_nr > -1) {
      _ekf._ekf.P(tropo_index(beacon_nr), tropo_index(beacon_nr)) = val;
    } else {
      const int start_idx = tropo_index(0);
      constexpr int step = 2;
      int idx = start_idx;
      for (int i = 0; i < _ekf.num_stations(); i++) {
        _ekf._ekf.P(idx, idx) = val;
        idx += step;
      }
    }
  }

  void set_rfoff_apriori(double val, int beacon_nr = -1) noexcept {
#ifdef DEBUG
    assert(beacon_nr >= -1 && beacon_nr < _ekf._num_stations);
#endif
    if (beacon_nr > -1) {
      _ekf._ekf.x(rfoff_index(beacon_nr)) = val;
    } else {
      const int start_idx = rfoff_index(0);
      constexpr int step = 2;
      int idx = start_idx;
      for (int i = 0; i < _ekf.num_stations(); i++) {
        _ekf._ekf.x(idx) = val;
        idx += step;
      }
    }
  }

  void set_rfoff_apriori_sigma(double val, int beacon_nr = -1) noexcept {
#ifdef DEBUG
    assert(beacon_nr >= -1 && beacon_nr < _ekf._num_stations);
#endif
    if (beacon_nr > -1) {
      _ekf._ekf.P(rfoff_index(beacon_nr), rfoff_index(beacon_nr)) = val;
    } else {
      const int start_idx = rfoff_index(0);
      constexpr int step = 2;
      int idx = start_idx;
      for (int i = 0; i < _ekf.num_stations(); i++) {
        _ekf._ekf.P(idx, idx) = val;
        idx += step;
      }
    }
  }

  void set_drag_coef_apriori(double val) noexcept {
    _ekf._ekf.x(drag_coef_index()) = val;
  }
  void set_drag_coef_apriori_sigma(double val) noexcept {
    _ekf._ekf.P(drag_coef_index(), drag_coef_index()) = val;
  }
};

} // namespace dso

#endif
