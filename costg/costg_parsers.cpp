#include "costg_utils.hpp"
#include <fstream>
#include <cstdio>
#include <charconv>
#include <cmath>
#include <cstring>

using costg::CostgAcc;
using costg::CostgExtState;
using costg::CostgQuat;

namespace {
  const char *skip_ws(const char *c) {
    while (*c && *c == ' ') ++c;
    return c;
  }
  int gravity_line(const char *line, CostgAcc &acc) {
    double d[4];
    const char *start = line;
    const char *end = line + std::strlen(line);
    for (int i=0; i<4; i++) {
      auto res = std::from_chars(skip_ws(start), end, d[i]);
      start = res.ptr;
      if (res.ec != std::errc{}) return 1;
    }
    double ip, fp;
    fp = std::modf(d[0], &ip);
    acc.gpst = dso::TwoPartDate(ip, fp);
    acc.ax = d[1];
    acc.ay = d[2];
    acc.az = d[3];
    return 0;
  }
  int estate_line(const char *line, CostgExtState &acc) {
    double d[10];
    const char *start = line;
    const char *end = line + std::strlen(line);
    for (int i=0; i<10; i++) {
      auto res = std::from_chars(skip_ws(start), end, d[i]);
      start = res.ptr;
      if (res.ec != std::errc{}) return 1;
    }
    double ip, fp;
    fp = std::modf(d[0], &ip);
    acc.gpst = dso::TwoPartDate(ip, fp);
    acc.x = d[1];
    acc.y = d[2];
    acc.z = d[3];
    acc.vx = d[4];
    acc.vy = d[5];
    acc.vz = d[6];
    acc.ax = d[7];
    acc.ay = d[8];
    acc.az = d[9];
    return 0;
  }
  int quaternion_line(const char *line, CostgQuat &acc) {
    double d[5];
    const char *start = line;
    const char *end = line + std::strlen(line);
    for (int i=0; i<5; i++) {
      auto res = std::from_chars(skip_ws(start), end, d[i]);
      start = res.ptr;
      if (res.ec != std::errc{}) return 1;
    }
    double ip, fp;
    fp = std::modf(d[0], &ip);
    acc.gpst = dso::TwoPartDate(ip, fp);
    acc.q = Eigen::Quaternion<double>(d[1],d[2],d[3],d[4]);
    return 0;
  }
}

int costg::parse_satellite_state(const char *fn, std::vector<CostgExtState> &acc) {
  acc.clear();

  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR. Failed to find input file %s\n", fn);
    return 1;
  }

  constexpr const int MC = 1024;
  char line[MC];

  while (fin.getline(line, MC)) {
    CostgExtState t;
    if (!estate_line(line, t)) acc.emplace_back(t);
  }

  return 0;
}

int costg::parse_gravity_field(const char *fn, std::vector<CostgAcc> &acc) {
  acc.clear();

  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR. Failed to find input file %s\n", fn);
    return 1;
  }

  constexpr const int MC = 512;
  char line[MC];

  while (fin.getline(line, MC)) {
    CostgAcc t;
    if (!gravity_line(line, t)) acc.emplace_back(t);
  }

  return 0;
}

int costg::parse_rotation_quaternions(const char *fn, std::vector<CostgQuat> &acc) {
  acc.clear();

  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR. Failed to find input file %s\n", fn);
    return 1;
  }

  constexpr const int MC = 512;
  char line[MC];

  while (fin.getline(line, MC)) {
    CostgQuat t;
    if (!quaternion_line(line, t)) acc.emplace_back(t);
  }

  return 0;
}
