#include "datetime/datetime_read.hpp"
#include "jason3_quaternions.hpp"
#include <charconv>
#include <cstdio>
#include <exception>
#include <fstream>
#include <stdexcept>

constexpr const int LINE_SZ = 256;

int resolve_jason3_body_quaternion_line(
    const char *line, dso::JasonBodyQuaternion &record) noexcept {
  // read in date (UTC)
  try {
    record.t = dso::strptime_ymd_hms<dso::nanoseconds>(line);
  } catch (std::exception &e) {
    fprintf(stderr,
            "[ERROR] Failed resolving date in quaternion file, lineL \'%s\' "
            "(traceback: %s)\n",
            line, __func__);
    return 1;
  }

  const char *last = line + 255;
  int error = 0;
  double qdata[4]; // temporary quaternion buffer

  // skip field UI1 and tab character
  const char *c = line+23+1;

  // UI1 -- 10 chars, unused value
  c += 10 + 1;

  // parse scalar part of quaternion
  auto pec = std::from_chars(c, last, qdata[0]);
  c = pec.ptr;
  if (*c++ != '\t' || pec.ec != std::errc{}) {
    ++error;
  }

  // parse vector part of quaternion
  for (int i = 0; i < 3; i++ && !error) {
    // UI2 -- 4 chars, unused value
    c += 4 + 1;
    // UI3 -- 10 chars, unused value
    c += 10 + 1;
    pec = std::from_chars(c, last, qdata[i + 1]);
    c = pec.ptr;
    if (*c++ != '\t' || pec.ec != std::errc{}) {
      ++error;
    }
  }

  // check for errors
  if (error) {
    fprintf(stderr,
            "[ERROR] Failed to resolve (body) Jason-3 quaternion line: \'%s\' "
            "(traceback: %s)\n",
            line, __func__);
    return 1;
  }

  // assign quaternion components
  // !! WARNING !!
  // it seems like the following line assignes first the vector part and then
  // the scalar:
  // record.quaternion = Eigen::Quaternion<double>(qdata);
  // this should do the trick:
  record.quaternion.w() = qdata[0];
  record.quaternion.x() = qdata[1];
  record.quaternion.y() = qdata[2];
  record.quaternion.z() = qdata[3];

  return 0;
}

dso::JasonQuaternionHunter::JasonQuaternionHunter(const char *body_fn)
    : bodyin(body_fn) {
  // read-in the first two quaternions
  if (bodyin.get_next(bodyq, 2)) {
    throw std::runtime_error(
        "ERROR JasonQuaternionHunter failed to caonstruct\n");
  }
}

int dso::JasonBodyQuaternionFile::get_next(
    dso::JasonBodyQuaternion &record) noexcept {
  char line[LINE_SZ];
  fin.getline(line, LINE_SZ);
  // skip commanet lines, if any
  while (line[0] == '#' && fin.good())
    fin.getline(line, LINE_SZ);
  return resolve_jason3_body_quaternion_line(line, record);
}

int dso::JasonBodyQuaternionFile::get_next(dso::JasonBodyQuaternion *records,
                                           int num_records) noexcept {
  char line[LINE_SZ];
  for (int i = 0; i < num_records; i++) {
    fin.getline(line, LINE_SZ);
    // skip comment lines, if any
    while (line[0] == '#' && fin.good())
      fin.getline(line, LINE_SZ);
    if (resolve_jason3_body_quaternion_line(line, records[i]))
      return i*10;
  }
  return 0;
}

int dso::JasonQuaternionHunter::set_at(
    const dso::datetime<dso::nanoseconds> &t) noexcept {
  // first, check if we are ok
  if (t >= bodyq[0].t && t < bodyq[1].t)
    return 0;

  // get next quaternion, hopefully we are ok now ...
  dso::JasonBodyQuaternion tmp;
  if (bodyin.get_next(tmp)) {
    fprintf(
        stderr,
        "[ERROR] Failed to get quaternion for requested date (traceback: %s)\n",
        __func__);
    return 1;
  }

  // check if we are ok with this pair ...
  if (t >= bodyq[1].t && t < tmp.t) {
    bodyq[0] = bodyq[1];
    bodyq[1] = tmp;

    return 0;
  } else {
    // else, find a suitable interval ....
    int error = 0;
    dso::JasonBodyQuaternion tmp2;
    do {
      error = bodyin.get_next(tmp2);
      bodyq[0] = tmp;
      bodyq[1] = tmp2;
      tmp = tmp2;
    } while (!error && !(t >= bodyq[0].t && t < bodyq[1].t));
    return error;
  }
}

//Eigen::Quaternion<double> slerp(const Eigen::Quaternion<double> &q1,
//                                const Eigen::Quaternion<double> &q2,
//                                double t) noexcept {
//  const double cOmega =
//      q1.w() * q2.w() + q1.x() * q2.x() + q1.y() * q2.y() + q1.z() * q2.z();
//  const double Omega = std::acos(cOmega);
//  const double s1mtO = std::sin((1e0 - t) * Omega);
//  const double stO = std::sin(t * Omega);
//  const double sO = std::sin(Omega);
//  Eigen::Quaternion<double> q = (q1 * s1mtO + q2 * stO) / sO;
//  return q.normalize();
//}

    int dso::JasonQuaternionHunter::get_at(
        const dso::datetime<dso::nanoseconds> &t,
        Eigen::Quaterniond &q) noexcept {
  if (set_at(t))
    return 1;

#ifdef DEBUG
  assert(t >= bodyq[0].t && t < bodyq[1].t);
#endif

  auto dts = bodyq[1].t.delta_sec(bodyq[0].t);
  const double dt_ab = dts.to_fractional_seconds(); // t2-t1
  dts = t.delta_sec(bodyq[0].t);
  const double dt_at = dts.to_fractional_seconds(); // t - t1
  const double dt = dt_at / dt_ab;
#ifdef DEBUG
  assert(dt >= 0e0 && dt <= 1e0);
#endif

  q = bodyq[0].quaternion.slerp(dt, bodyq[1].quaternion);
  q.normalize();
  // q = slerp(bodyq[0].quaternion, bodyq[1].quaternion, dt);
  return 0;
}
