#ifndef __FILTER_IMPLEMENTATION_MODELS_HPP__
#define __FILTER_IMPLEMENTATION_MODELS_HPP__

#include "eigen3/Eigen/Eigen"
#include "filters/ekf.hpp"
#include "datetime/dtcalendar.hpp"
#ifdef DEBUG
#include <cstdio>
#include <cassert>
#endif

namespace dso {

enum class BeaconClockModel : char {None, Linear};

template<int Np, BeaconClockModel bcm>
struct EkfFilterImpl {
    ExtendedKalmanFilter<dso::nanoseconds> _ekf;
    dso::datetime<dso::nanoseconds> _ref_epoch;
    int _num_stations{0};

    static std::size_t num_params(int num_stations) noexcept {
        return 
        6    // satellite state vector
        + Np // satellite force model parameters
        + num_stations // zenith tropospheric wet delay, per beacon
        // number of relative clock offset parameters, per beacon
        + num_stations * ((bcm==BeaconClockModel::None) ? (1) : (2));
    }

    int num_params() const noexcept {return EkfFilterImpl<Np,bcm>::num_params(_num_stations);}

    int num_stations() const noexcept {return _num_stations;}

    Eigen::MatrixXd estimates() const noexcept {
        return _ekf.x;
    }

    Eigen::MatrixXd &P() noexcept {
        return _ekf.P;
    }

    EkfFilterImpl(int num_stations,
              const dso::datetime<dso::nanoseconds> ref_epoch) noexcept
        : _ekf(EkfFilterImpl<Np,bcm>::num_params(num_stations)), _ref_epoch(ref_epoch),
          _num_stations(num_stations) {}
    
    void time_update(const dso::datetime<dso::nanoseconds> &tk,
                     const Eigen::VectorXd &xk,
                     const Eigen::MatrixXd &phi) noexcept {
      _ekf.time_update(tk, xk, phi);
    }

    void observation_update(double z, double g, double sigma,
                            const Eigen::VectorXd &H) noexcept {
      _ekf.observation_update(z, g, sigma, H);
    }

}; // ekfFilterImpl

template<int Np, BeaconClockModel bcm>
struct EkfFilter {};

template<int Np>
struct EkfFilter<Np, BeaconClockModel::None> {
    EkfFilterImpl<Np, BeaconClockModel::None> _ekf;

    EkfFilter(int num_stations,
              const dso::datetime<dso::nanoseconds> ref_epoch) noexcept
        : _ekf(num_stations, ref_epoch) {}
    
    std::size_t num_params() const noexcept {
        return _ekf.num_params();
    }

    Eigen::MatrixXd estimates() const noexcept {
        return _ekf.estimates();
    }
    Eigen::MatrixXd &P() noexcept {
        return _ekf.P();
    }

    double tropo_estimate(int beacon_nr) const noexcept {
        return _ekf._ekf.x(6 + Np + beacon_nr * 2);
    }
    double &tropo_estimate(int beacon_nr) noexcept {
        return _ekf._ekf.x(6 + Np + beacon_nr * 2);
    }
    double rfoff_estimate(int beacon_nr) const noexcept {
        return _ekf._ekf.x(6 + Np + beacon_nr * 2 + 1);
    }
    double &rfoff_estimate(int beacon_nr) noexcept {
        return _ekf._ekf.x(6 + Np + beacon_nr * 2 + 1);
    }

    void time_update(const dso::datetime<dso::nanoseconds> &tk,
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
      constexpr int start_idx = 6 + Np;
      if (beacon_nr > -1) {
        _ekf._ekf.x(start_idx + beacon_nr*2) = val;
      } else {
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
      constexpr int start_idx = 6 + Np;
      if (beacon_nr > -1) {
        _ekf._ekf.P(start_idx + beacon_nr*2, start_idx + beacon_nr*2) = val;
      } else {
        constexpr int step = 2;
        int idx = start_idx;
        for (int i = 0; i < _ekf.num_stations(); i++) {
          _ekf._ekf.P(idx,idx) = val;
          idx += step;
        }
      }
    }

    void set_rfoff_apriori(double val, int beacon_nr = -1) noexcept {
#ifdef DEBUG
      assert(beacon_nr >= -1 && beacon_nr < _ekf._num_stations);
#endif
      constexpr int start_idx = 6 + Np + 1;
      if (beacon_nr > -1) {
        _ekf._ekf.x(start_idx + beacon_nr*2) = val;
      } else {
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
      constexpr int start_idx = 6 + Np + 1;
      if (beacon_nr > -1) {
        _ekf._ekf.P(start_idx + beacon_nr*2, start_idx + beacon_nr*2) = val;
      } else {
        constexpr int step = 2;
        int idx = start_idx;
        for (int i = 0; i < _ekf.num_stations(); i++) {
          _ekf._ekf.P(idx,idx) = val;
          idx += step;
        }
      }
    }

};

}// dso

#endif