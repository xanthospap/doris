#include "datetime/utcdates.hpp"
#include "jason3_quaternions.hpp"
#include <charconv>
#include <cstdio>
#include <datetime/dtcalendar.hpp>
#include <datetime/dtfund.hpp>
#include <exception>
#include <fstream>
#include <stdexcept>

namespace {
constexpr const int LINE_SZ = 512;

int resolve_jason3_body_quaternion_line(
    const char *line, dso::JasonBodyQuaternion &record) noexcept {
  // read in date (UTC) and transform to TAI
  try {
    auto utc = dso::utc_strptime_ymd_hms(line);
    record.tai_mjd = dso::utc2tai(utc);
    // printf("Note that we resolved ITC date %s to %.2f + %.15f\n", line, record.tai_mjd._big, record.tai_mjd._small);
  } catch (std::exception &e) {
    fprintf(stderr,
            "[ERROR] Failed resolving date in quaternion file, line \'%s\' "
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
  if (!fin.getline(line, LINE_SZ)) {
    fprintf(stderr, "[ERROR] Failed to read line from quaternion file! (traceback: %s)\n", __func__);
    return 2;
  }
  
  // skip commanet lines, if any
  while (line[0] == '#' && fin.good())
    fin.getline(line, LINE_SZ);
  
  return resolve_jason3_body_quaternion_line(line, record);
}

int dso::JasonBodyQuaternionFile::get_next(dso::JasonBodyQuaternion *records,
                                           int num_records) noexcept {
  char line[LINE_SZ];
  for (int i = 0; i < num_records; i++) {
    if (!fin.getline(line, LINE_SZ)) {
      fprintf(stderr, "[ERROR] Failed to read line from quaternion file! (traceback: %s)\n", __func__);
      return 2;
    }

    // skip comment lines, if any
    while (line[0] == '#' && fin.good())
      fin.getline(line, LINE_SZ);
    
    if (resolve_jason3_body_quaternion_line(line, records[i]))
      return i*10;
  }

  return 0;
}

int dso::JasonQuaternionHunter::set_at(const dso::TwoPartDate &tai_mjd, int &index) noexcept {
  // first, check if we are ok, i.e. the time requested is within the buffered
  // interval
  index = this->find_interval(tai_mjd);
  // printf("Index returned: %d\n", index);
  assert(index >= 0);
  if (index < NumQuaternionsInBuffer) {
    return 0;
  }

#ifdef DEBUG
  assert(tmp.tai_mjd > bodyq[NumQuaternionsInBuffer - 1].tai_mjd);
#endif

  // get next quaternion, hopefully we are ok now ...
  // new quternion temporarilly stored in tmp
  dso::JasonBodyQuaternion tmp;
  if (bodyin.get_next(tmp)) {
    fprintf(
        stderr,
        "[ERROR] Failed to get quaternion for requested date (traceback: %s)\n",
        __func__);
    return 1;
  }

  // check if we are ok with this pair ...
  if ((tai_mjd >= bodyq[NumQuaternionsInBuffer - 1].tai_mjd) &&
      (tai_mjd < tmp.tai_mjd)) {
    //  left shift quaternions and add the new in the last index
    left_shift();
    bodyq[NumQuaternionsInBuffer - 1] = tmp;
    index = NumQuaternionsInBuffer - 2;

    return 0;
  } else {
    if (tai_mjd < bodyq[0].tai_mjd) {
      fprintf(stderr,
              "ERROR Requested quaternion for an epoch prior to current "
              "interval, t=%.15f\n",
              tai_mjd.mjd());
      fprintf(stderr, "      Current interval spans: %.15f to %.15f\n",
              bodyq[0].tai_mjd.mjd(), bodyq[1].tai_mjd.mjd());
      return 1;
    }
    //printf("Hunting failed, next quaternions was at: %.1f + %.15f\n", tmp.tai_mjd._big, tmp.tai_mjd._small);
    //printf("Rejected beacuse: %1d %1d\n",(tai_mjd >= bodyq[NumQuaternionsInBuffer - 1].tai_mjd), (tai_mjd < tmp.tai_mjd));
    //printf("Again, that is: %.1f + %.15f < %.1f + %.15f\n", tai_mjd._big, tai_mjd._small, tmp.tai_mjd._big, tmp.tai_mjd._small);
    // else, find a suitable interval ....
    int error = 0;
    left_shift();
    bodyq[NumQuaternionsInBuffer - 1] = tmp;

    do {
      error = bodyin.get_next(tmp);
      left_shift();
      bodyq[NumQuaternionsInBuffer - 1] = tmp;
      // printf("\tnext quaternions is at: %.1f + %.15f\n", tmp.tai_mjd._big, tmp.tai_mjd._small);
      // bodyq[0] = tmp;
      // bodyq[1] = tmp2;
      // tmp = tmp2;
    } while (!error && !((tai_mjd >= bodyq[NumQuaternionsInBuffer - 2].tai_mjd) &&
                         (tai_mjd < bodyq[NumQuaternionsInBuffer - 1].tai_mjd)));
    if (error) {
      fprintf(stderr,
              "[ERROR] Failed hunting quaternion for epoch: %.15f. Last found "
              "was for %.15f. error=%d (traceback: %s)\n",
              tai_mjd.mjd(), bodyq[NumQuaternionsInBuffer - 1].tai_mjd.mjd(),
              error, __func__);
    }

    return (error) ? (error) : ((index = NumQuaternionsInBuffer - 2) == -1);
  }
}

int dso::JasonQuaternionHunter::get_at(const dso::TwoPartDate &tai_mjd,
                                       Eigen::Quaterniond &q) noexcept {
  constexpr const dso::DateTimeDifferenceType FD = dso::DateTimeDifferenceType::FractionalDays;
  int index;
  if (set_at(tai_mjd, index))
    return 99;

#ifdef DEBUG
  assert(index>=0 && index < NumQuaternionsInBuffer-1);
  assert(tai_mjd >= bodyq[index].tai_mjd && tai_mjd < bodyq[index+1].tai_mjd);
#endif

  const int i0 = index;
  const int i1 = index+1;
  
  const double dt_ab = bodyq[i1].tai_mjd.diff<FD>(bodyq[i0].tai_mjd);
  const double dt_at = tai_mjd.diff<FD>(bodyq[i0].tai_mjd);
  const double dt = dt_at / dt_ab;
#ifdef DEBUG
  assert(dt >= 0e0 && dt <= 1e0);
#endif
  // Eigen way
  q = bodyq[i0].quaternion.slerp(dt, bodyq[i1].quaternion);
  q.normalize();
  // q = slerp(bodyq[0].quaternion, bodyq[1].quaternion, dt);
  return 0;
}
