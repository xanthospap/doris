#include "datetime/utcdates.hpp"
#include "jason3_quaternions.hpp"
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

/* Resolve a solar array orientation line. Note that the time record is
 * tranlated to TAI (from UTC)
 */
int resolve_jason3_solar_array_line(
    const char *line, dso::JasonSolarArrayRotation &record) noexcept {
  /* read in date (UTC) and transform to TAI */
  try {
    auto utc = dso::utc_strptime_ymd_hms(line);
    record.tai_mjd = utc.utc2tai();
  } catch (std::exception &e) {
    fprintf(stderr,
            "[ERROR] Failed resolving date in solar array file, line \'%s\' "
            "(traceback: %s)\n",
            line, __func__);
    return 1;
  }

  const char *last = line + std::strlen(line);
  int error = 0;
  double angle; /* temporary angle buffer */

  /* skip field UI1 and space/tab character */
  const char *c = line + 23 + 1;
  /* UI1 -- 10 chars, unused value */
  c += 10 + 1;

  auto pec = std::from_chars(skip_ws(c), last, angle);
  c = pec.ptr;
  if (!std::isspace(*c) || pec.ec != std::errc{}) {
    ++error;
  }
  record.left_array = angle;

  /* Useless parameter, format (I10) */
  c += 10;
  pec = std::from_chars(skip_ws(c), last, angle);
  if (!std::isspace(*c) || pec.ec != std::errc{}) {
    ++error;
  }
  record.right_array = angle;

  /* check for errors */
  if (error) {
    fprintf(stderr,
            "[ERROR] Failed to resolve Jason-3 solar array line: \'%s\' "
            "(traceback: %s)\n",
            line, __func__);
    return 1;
  }

  return 0;
}
} /* unnamed namespace */

dso::JasonSolarArrayHunter::JasonSolarArrayHunter(const char *body_fn)
    : bodyin(body_fn) {
  /* read-in the first NumQuaternionsInBuffer quaternions */
  if (bodyin.get_next(rots, NumQuaternionsInBuffer)) {
    throw std::runtime_error(
        "ERROR JasonSolarArrayHunter failed to construct\n");
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
int dso::JasonSolarArrayHunter::find_interval(
    const dso::TwoPartDate &tai_mjd) const noexcept {

  /* start searching from the top, aka from last element */
  int qindex = NumQuaternionsInBuffer - 2;
  for (int i = qindex; i >= 0; --i) {
    if ((tai_mjd >= rots[i].tai_mjd) && (tai_mjd < rots[i + 1].tai_mjd)) {
      return i;
    }
  }

  /* tai_mjd is out of bounds, prior to first record in buffer */
  if (tai_mjd < rots[0].tai_mjd) {
    fprintf(stderr,
            "[WRNNG] First rotations in buffer is at %.9f, requested epoch "
            "%.9f (traceback: %s)\n",
            rots[0].tai_mjd.as_mjd(), tai_mjd.as_mjd(), __func__);
    return std::numeric_limits<int>::min();
  }

  /* tai_mjd is out of bounds, after the last record in buffer */
  if (tai_mjd >= rots[NumQuaternionsInBuffer - 1].tai_mjd)
    return std::numeric_limits<int>::max();

  /* we should never reach this point */
  fprintf(stderr,
          "[ERROR] Hit wall while searching for solar array angles: requested date: "
          "%.9f, buffered: %.9f to %.9f (traceback: %s)\n",
          tai_mjd.as_mjd(), rots[0].tai_mjd.as_mjd(),
          rots[NumQuaternionsInBuffer - 1].tai_mjd.as_mjd(), __func__);
  assert(1 == 2);
}

/* Read in the next record off from the Solar Array input stream and store
 * it in the passed in record instance
 */
dso::iStatus dso::JasonSolarArrayFile::get_next(
    dso::JasonSolarArrayRotation &record) noexcept {
  char line[LINE_SZ];
  /* read next line */
  if (!fin.getline(line, LINE_SZ)) {
    fprintf(
        stderr,
        "[ERROR] Failed to read line from solar array file! (traceback: %s)\n",
        __func__);
    return dso::iStatus(2);
  }
  /* skip comment lines, if any */
  while (line[0] == '#' && fin.good())
    fin.getline(line, LINE_SZ);
  /* resolve line into record */
  return resolve_jason3_solar_array_line(line, record);
}

/* Read in the next num_records records off from the Solar Array input
 * stream and store them in the passed-in record instance in indexes
 * [0,num_records)
 */
dso::iStatus
dso::JasonSolarArrayFile::get_next(dso::JasonSolarArrayRotation *records,
                                       int num_records) noexcept {
  char line[LINE_SZ];
  int error = 0;
  int records_extracted = 0;

  /* iteratively extract num_records records ... */
  while (records_extracted < num_records && !error) {
    if (!fin.getline(line, LINE_SZ)) {
      fprintf(
          stderr,
          "[ERROR] Failed to read line from solar array file! (traceback: %s)\n",
          __func__);
      return dso::iStatus(2);
    }

    /* skip comment lines, if any */
    if (line[0] != '#') {
      error += (!fin.good());
      error +=
          (resolve_jason3_solar_array_line(line, records[records_extracted]));
      records_extracted += (!error) * 1;
    }
  } /* num_records records extracted */

  return dso::iStatus(error);
}

dso::iStatus dso::JasonSolarArrayHunter::set_at(const dso::TwoPartDate &tai_mjd,
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
            "[WRNNG] Rewinding solar array stream to find suitable interval "
            "(traceback: %s)\n",
            __func__);
    /* rewind stream */
    bodyin.fin.clear();
    bodyin.fin.seekg(0, std::ios::beg);
    /* collect first NumQuaternionsInBuffer angles */
    if (bodyin.get_next(rots, NumQuaternionsInBuffer)) {
      fprintf(stderr,
              "[ERROR] Failed to collect initial solar array angles after rewiding! "
              "(traceback: %s)\n",
              __func__);
      return dso::iStatus(99);
    }
    /* call function again, from the top of the file ... */
    return set_at(tai_mjd, index);
  }

  /* try getting next quaternion, hopefully we are ok now ... */
  dso::JasonSolarArrayRotation tmp;
  if (bodyin.get_next(tmp)) {
    fprintf(
        stderr,
        "[ERROR] Failed to get solar array angles for requested date (traceback: %s)\n",
        __func__);
    return dso::iStatus(1);
  }
#ifdef DEBUG
  assert(tmp.tai_mjd > rots[NumQuaternionsInBuffer - 1].tai_mjd);
#endif
  /* push back to buffer */
  left_shift();
  rots[NumQuaternionsInBuffer - 1] = tmp;
  /* check if we are ok with this pair ... */
  if ((tai_mjd >= rots[NumQuaternionsInBuffer - 2].tai_mjd) &&
      (tai_mjd < rots[NumQuaternionsInBuffer - 1].tai_mjd)) {
    index = NumQuaternionsInBuffer - 2;
    return dso::iStatus(0);
  }

  /* nope, need to read in untill we match the time interval */
  return this->read_untill_buffered(tai_mjd, index);
}

dso::iStatus dso::JasonSolarArrayHunter::read_untill_buffered(
    const dso::TwoPartDate &tai_mjd, int &index) noexcept {
  constexpr const int N = 2;
  dso::iStatus error = dso::iStatus(0);
  do {
    left_shift(N);
    error = bodyin.get_next(rots + NumQuaternionsInBuffer - N, N);
    index = this->find_interval(tai_mjd);
  } while (!error && (index != std::numeric_limits<int>::min() &&
                      index != std::numeric_limits<int>::max()));

  return error;
}

dso::iStatus dso::JasonSolarArrayHunter::get_at(const dso::TwoPartDate &tai_mjd,
                                                double &left_angle,
                                                double &right_angle) noexcept {
  constexpr const dso::DateTimeDifferenceType FD =
      dso::DateTimeDifferenceType::FractionalDays;
  int index;

  /* find the record for the given epoch */
  if (set_at(tai_mjd, index))
    return dso::iStatus(99);

#ifdef DEBUG
  assert(index >= 0 && index < NumQuaternionsInBuffer - 1);
  assert(tai_mjd >= rots[index].tai_mjd && tai_mjd < rots[index + 1].tai_mjd);
#endif

  const int i0 = index;     /* t1 */
  const int i1 = index + 1; /* t2 */

  const double dt_ab = rots[i1].tai_mjd.diff<FD>(
      rots[i0].tai_mjd); /* t2-t1 as fractional days */
  const double dt_at =
      tai_mjd.diff<FD>(rots[i0].tai_mjd); /* t-t1 as fractional days */
  const double dt_bt = rots[i1].tai_mjd.diff<FD>(tai_mjd); /* t2-t */

  /* linear interpolation */
  left_angle = rots[i0].left_array * (dt_bt / dt_ab) +
               rots[i1].left_array * (dt_at / dt_ab);
  right_angle = rots[i0].right_array * (dt_bt / dt_ab) +
                rots[i1].right_array * (dt_at / dt_ab);

  /* all done */
  return dso::iStatus(0);
}
