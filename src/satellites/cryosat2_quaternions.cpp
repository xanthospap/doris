#include "datetime/datetime_read.hpp"
#include "jason3_quaternions.hpp"
#include <charconv>
#include <cstdio>
#include <exception>
#include <fstream>
#include <stdexcept>

constexpr const int LINE_SZ = 256;

int resolve_cryosat2_body_quaternion_line(
    const char *line, dso::JasonBodyQuaternion &record) noexcept {
  // read in date (UTC)
  try {
    record.t = dso::strptime_ymd_hms<dso::nanoseconds>(line);
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

  // skip date (not time)
  const char *c = line+11;
  // first whitespace after time
  while (*c++ != ' ') ; 

  // parse vector part of quaternion
  for (int i = 0; i < 3; i++ && !error) {
    auto pec = std::from_chars(c, last, qdata[i + 1]);
    c = pec.ptr;
    if (*c++ != ' ' || pec.ec != std::errc{}) {
      ++error;
    } 
    while (*c++ == ' ') ; // skip whitespace characters ... 
  }
  
  // parse scalar part of quaternion
  auto pec = std::from_chars(c, last, qdata[0]);
  c = pec.ptr;
  if (*c++ != ' ' || pec.ec != std::errc{}) {
    ++error;
  }

  // check for errors
  if (error) {
    fprintf(stderr,
            "[ERROR] Failed to resolve (body) Cryosat-2 quaternion line: \'%s\' "
            "(traceback: %s)\n",
            line, __func__);
    return 1;
  }

  // assign quaternion components
  record.quaternion.w() = qdata[0];
  record.quaternion.x() = qdata[1];
  record.quaternion.y() = qdata[2];
  record.quaternion.z() = qdata[3];

  return 0;
}