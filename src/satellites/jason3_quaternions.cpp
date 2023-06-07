#include "jason3_quaternions.hpp"
#include "datetime/utcdates.hpp"
#include <cctype>
#include <charconv>
#include <cstdio>
#include <exception>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace {
constexpr const int LINE_SZ = 512;

const char *skip_ws(const char *str) noexcept {
  while (*str && std::isspace(*str))
    ++str;
  return str;
}

/* Resolve a body-quaternion line. Note that the time record is tranlated
 * to TAI (from UTC)
 */
int resolve_jason3_body_quaternion_line(
    const char *line, dso::JasonBodyQuaternion &record) noexcept {
  /* read in date (UTC) and transform to TAI */
  try {
    auto utc = dso::utc_strptime_ymd_hms(line);
    record.tai_mjd = utc.utc2tai();
  } catch (std::exception &e) {
    fprintf(stderr,
            "[ERROR] Failed resolving date in quaternion file, line \'%s\' "
            "(traceback: %s)\n",
            line, __func__);
    return 1;
  }

  const char *last = line + std::strlen(line);
  int error = 0;
  double qdata[4]; /* temporary quaternion buffer */

  /* skip field UI1 and space/tab character */
  const char *c = line + 23 + 1;
  /* UI1 -- 10 chars, unused value */
  c += 10 + 1;

  /* parse scalar part of quaternion */
  auto pec = std::from_chars(skip_ws(c), last, qdata[0]);
  c = pec.ptr;
  if (!std::isspace(*c) || pec.ec != std::errc{}) {
    fprintf(stderr, "[ERROR] Failed parsing scalar part of quaternion line!\n");
    ++error;
  }

  /* parse vector part of quaternion */
  for (int i = 0; i < 3 && (!error); i++) {
    /* 4 chars + 1 tab + 10 chars */
    c += 15;
    pec = std::from_chars(skip_ws(c), last, qdata[i + 1]);
    c = pec.ptr;
    if (!std::isspace(*c) || pec.ec != std::errc{}) {
      fprintf(stderr,
              "[ERROR] Failed parsing imaginary part #%d of quaternion line!\n",
              i);
      ++error;
    }
  }

  /* check for errors */
  if (error) {
    fprintf(stderr,
            "[ERROR] Failed to resolve (body) Jason-3 quaternion line: \'%s\' "
            "(traceback: %s)\n",
            line, __func__);
    return 1;
  }

  /* assign quaternion components
   * !! WARNING !!
   * it seems like the following line assignes first the vector part and then
   * the scalar:
   * record.quaternion = Eigen::Quaternion<double>(qdata);
   * this should do the trick:
   */
  record.quaternion.w() = qdata[0];
  record.quaternion.x() = qdata[1];
  record.quaternion.y() = qdata[2];
  record.quaternion.z() = qdata[3];

  return 0;
}
} /* unnamed namespace */

dso::JasonQuaternionHunter::JasonQuaternionHunter(const char *body_fn)
    : bodyin(body_fn) {
  /* read-in the first NumQuaternionsInBuffer quaternions */
  if (bodyin.get_next(bodyq, NumQuaternionsInBuffer)) {
    throw std::runtime_error(
        "ERROR JasonQuaternionHunter failed to construct\n");
  }
}

/* Find an interval within the buffered quaternions for which :
 * t >= q[i] && t < q[i+1]
 * (where q is the buffered quaternion array) and return i.
 * If no such interval is found and t < q[0], then
 * std::numeric_limits<int>::min() is returned. If no such interval is found
 * and t >= q[maxNumBufferedQuaternons], then std::numeric_limits<int>::max()
 * is returned
 */
int dso::JasonQuaternionHunter::find_interval(
    const dso::TwoPartDate &tai_mjd) const noexcept {

  /* start searching from the top, aka from last element */
  int qindex = NumQuaternionsInBuffer - 2;
  for (int i = qindex; i >= 0; --i) {
    if ((tai_mjd >= bodyq[i].tai_mjd) && (tai_mjd < bodyq[i + 1].tai_mjd)) {
      return i;
    }
  }

  /* tai_mjd is out of bounds, prior to first record in buffer */
  if (tai_mjd < bodyq[0].tai_mjd) {
    fprintf(stderr,
            "[WRNNG] First quaternion in buffer is at %.9f, requested epoch "
            "%.9f (traceback: %s)\n",
            bodyq[0].tai_mjd.as_mjd(), tai_mjd.as_mjd(), __func__);
    return std::numeric_limits<int>::min();
  }

  /* tai_mjd is out of bounds, after the last record in buffer */
  if (tai_mjd >= bodyq[NumQuaternionsInBuffer - 1].tai_mjd)
    return std::numeric_limits<int>::max();

  /* we should never reach this point */
  fprintf(stderr,
          "[ERROR] Hit wall while searching for quaternion: requested date: "
          "%.9f, buffered: %.9f to %.9f (traceback: %s)\n",
          tai_mjd.as_mjd(), bodyq[0].tai_mjd.as_mjd(),
          bodyq[NumQuaternionsInBuffer - 1].tai_mjd.as_mjd(), __func__);
  assert(1 == 2);
}

/* Read in the next record off from the Body Quaternion input stream and store
 * it in the passed in record instance
 */
dso::iStatus dso::JasonBodyQuaternionFile::get_next(
    dso::JasonBodyQuaternion &record) noexcept {
  char line[LINE_SZ];
  /* read next line */
  if (!fin.getline(line, LINE_SZ)) {
    fprintf(
        stderr,
        "[ERROR] Failed to read line from quaternion file! (traceback: %s)\n",
        __func__);
    return dso::iStatus(2);
  }
  /* skip comment lines, if any */
  while (line[0] == '#' && fin.good())
    fin.getline(line, LINE_SZ);
  /* resolve line into record */
  return resolve_jason3_body_quaternion_line(line, record);
}

/* Read in the next num_records records off from the Body Quaternion input
 * stream and store them in the passed in record instance in indexes
 * [0,num_records)
 */
dso::iStatus
dso::JasonBodyQuaternionFile::get_next(dso::JasonBodyQuaternion *records,
                                       int num_records) noexcept {
  char line[LINE_SZ];
  int error = 0;
  int records_extracted = 0;

  /* iteratively extract num_records records ... */
  while (records_extracted < num_records && !error) {
    if (!fin.getline(line, LINE_SZ)) {
      fprintf(
          stderr,
          "[ERROR] Failed to read line from quaternion file! (traceback: %s)\n",
          __func__);
      return dso::iStatus(2);
    }

    /* skip comment lines, if any */
    if (line[0] != '#') {
      error += (!fin.good());
      error += (resolve_jason3_body_quaternion_line(
          line, records[records_extracted]));
      records_extracted += (!error) * 1;
    }
  } /* num_records records extracted */

  return dso::iStatus(error);
}

dso::iStatus dso::JasonQuaternionHunter::set_at(const dso::TwoPartDate &tai_mjd,
                                                int &index) noexcept {
  /* first, check if we are ok, i.e. the time requested is within the
   * buffered interval
   */
  index = this->find_interval(tai_mjd);
  if (index != std::numeric_limits<int>::min() &&
      index != std::numeric_limits<int>::max()) {
    return dso::iStatus(0);
  }

  /* Must restart searching from the top! */
  if (index == std::numeric_limits<int>::min()) {
    fprintf(stderr,
            "[WRNNG] Rewinding quaternion stream to find suitable interval "
            "(traceback: %s)\n",
            __func__);
    /* rewind stream */
    bodyin.fin.clear();
    bodyin.fin.seekg(0, std::ios::beg);
    /* collect first NumQuaternionsInBuffer quaternions */
    if (bodyin.get_next(bodyq, NumQuaternionsInBuffer)) {
      fprintf(stderr,
              "[ERROR] Failed to collect initial quaternions after rewiding! "
              "(traceback: %s)\n",
              __func__);
      return dso::iStatus(99);
    }
    /* call function again, from the top of the file ... */
    return set_at(tai_mjd, index);
  }

  /* try getting next quaternion, hopefully we are ok now ... */
  dso::JasonBodyQuaternion tmp;
  if (bodyin.get_next(tmp)) {
    fprintf(
        stderr,
        "[ERROR] Failed to get quaternion for requested date (traceback: %s)\n",
        __func__);
    return dso::iStatus(1);
  }
#ifdef DEBUG
  assert(tmp.tai_mjd > bodyq[NumQuaternionsInBuffer - 1].tai_mjd);
#endif
  /* push back to buffer */
  left_shift();
  bodyq[NumQuaternionsInBuffer - 1] = tmp;
  /* check if we are ok with this pair ... */
  if ((tai_mjd >= bodyq[NumQuaternionsInBuffer - 2].tai_mjd) &&
      (tai_mjd < bodyq[NumQuaternionsInBuffer - 1].tai_mjd)) {
    index = NumQuaternionsInBuffer - 2;
    return dso::iStatus(0);
  }

  /* nope, need to read in untill we match the time interval */
  return this->read_untill_buffered(tai_mjd, index);
}

dso::iStatus dso::JasonQuaternionHunter::read_untill_buffered(
    const dso::TwoPartDate &tai_mjd, int &index) noexcept {
  constexpr const int N = 2;
  dso::iStatus error = dso::iStatus(0);
  do {
    left_shift(N);
    error = bodyin.get_next(bodyq + NumQuaternionsInBuffer - N, N);
    index = this->find_interval(tai_mjd);
  } while (!error && index == std::numeric_limits<int>::max());

  return error;
}

dso::iStatus
dso::JasonQuaternionHunter::get_at(const dso::TwoPartDate &tai_mjd,
                                   Eigen::Quaterniond &q) noexcept {
  constexpr const dso::DateTimeDifferenceType FD =
      dso::DateTimeDifferenceType::FractionalDays;
  int index;

  /* find the quaternion record for the given epoch */
  if (set_at(tai_mjd, index))
    return dso::iStatus(99);

#ifdef DEBUG
  if (index<0 || index >= NumQuaternionsInBuffer - 1) {
    printf("Buffered quaternions:\n");
    for (int i=0; i<NumQuaternionsInBuffer; i++) {
      printf("[%2d] %.9f\n", i, bodyq[i].tai_mjd.as_mjd());
    }
    printf("Requested quaternion at %.9f, index found = %d\n", tai_mjd.as_mjd(), index);
  }
  assert(index >= 0 && index < NumQuaternionsInBuffer - 1);
  assert(tai_mjd >= bodyq[index].tai_mjd && tai_mjd < bodyq[index + 1].tai_mjd);
#endif

  const int i0 = index;     /* t1 */
  const int i1 = index + 1; /* t2 */

  const double dt_ab = bodyq[i1].tai_mjd.diff<FD>(
      bodyq[i0].tai_mjd); /* t2-t as fractional days */
  const double dt_at =
      tai_mjd.diff<FD>(bodyq[i0].tai_mjd); /* t-t1 as fractional days */
  const double dt = dt_at / dt_ab;
#ifdef DEBUG
  assert(dt >= 0e0 && dt <= 1e0);
#endif

  /* SLERP interpolation */
  q = bodyq[i0].quaternion.slerp(dt, bodyq[i1].quaternion);

  /* normalize the quaternion */
  q.normalize();

  /* all done */
  return dso::iStatus(0);
}
