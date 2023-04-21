#ifndef __COSTG_HELPER_FUNCTIONS_HPP__
#define __COSTG_HELPER_FUNCTIONS_HPP__

#include "datetime/dtcalendar.hpp"
#include "eigen3/Eigen/Eigen"
#include "eigen3/Eigen/Geometry"
#include <vector>

namespace costg {
struct CostgAcc {
  /* gps time */
  dso::TwoPartDate gpst;
  /* acceleration components in [m/s^2] */
  double ax, ay, az;

  Eigen::Matrix<double, 3, 1> a() {
    return Eigen::Matrix<double, 3, 1>({ax, ay, az});
  }
}; /* CostgAcc */

struct CostgExtState {
  /* gps time */
  dso::TwoPartDate gpst;
  /* position components in [m] */
  double x, y, z;
  Eigen::Matrix<double, 3, 1> r() {
    return Eigen::Matrix<double, 3, 1>({x, y, z});
  }
  /* velocity components in [m/sec] */
  double vx, vy, vz;
  Eigen::Matrix<double, 3, 1> v() {
    return Eigen::Matrix<double, 3, 1>({vx, vy, vz});
  }
  /* acceleration components in [m/s^2] */
  double ax, ay, az;
  Eigen::Matrix<double, 3, 1> a() {
    return Eigen::Matrix<double, 3, 1>({ax, ay, az});
  }
}; /* CostgExtState */

struct CostgQuat {
  /* gps time */
  dso::TwoPartDate gpst;
  /* position components in [m] */
  Eigen::Quaternion<double> q;
};

int parse_gravity_field(const char *fn, std::vector<CostgAcc> &acc);
int parse_satellite_state(const char *fn, std::vector<CostgExtState> &acc);
int parse_rotation_quaternions(const char *fn, std::vector<CostgQuat> &quats);

inline dso::TwoPartDate gps2tai(const dso::TwoPartDate &gps) {
  return dso::TwoPartDate(gps.big(), gps.small() + 19e0 / dso::sec_per_day);
}

} /* namespace costg */
#endif
