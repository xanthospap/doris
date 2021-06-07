#include "sp3.hpp"
#include <stdexcept>
#include <cstdio>
#ifdef DEBUG
#include "ggdatetime/datetime_write.hpp"
#include <iostream>
#endif

/// Max record characters (for a navigation data block)
constexpr int MAX_RECORD_CHARS{128};

/// Bad or absent clock values areto be set to 999999.999999.  The six integer
/// nines are required, whereasthe fractional part nines are optional.
constexpr double SP3_MISSING_CLK_VALUE{999999.e0};

/// @brief Check if given substring is empty (aka whitespace only).
/// Substring to check is str with indexes [0,count)
/// @param[in] str Start of string
/// @param[in] count Number of chars to consider
/// @return true id substring is whitespace only, false otherwise
bool substr_is_empty(const char *str, std::size_t count) noexcept {
  std::size_t idx = 0;
  while (idx < count && str[idx] == ' ')
    ++idx;
  return idx == count;
}

/// @brief Resolve an Epoch Header Record line
/// @param[in] line An Epoch Header Record to be resolved
/// @param[out] t The epoch resolved from the input line
/// @return Anything other than 0 denotes an error
int ids::Sp3c::resolve_epoch_line(
    ngpt::datetime<ngpt::microseconds> &t) noexcept {
  char line[MAX_RECORD_CHARS];
  char *end;
  const char *start;

  __istream.getline(line, MAX_RECORD_CHARS);
  if (line[0] != '*' || line[1] != ' ')
    return 2;

  int date[5];
  date[0] = std::strtol(line + 3, &end, 10);
  if (!date[0] || errno == ERANGE) {
    errno = 0;
    return 5;
  }
  start = line + 8;
  for (int i = 1; i < 5; i++) {
    date[i] = std::strtol(start, &end, 10);
    if (errno == ERANGE || end == start) {
      errno = 0;
      return 5 + i;
    }
    start += 3;
  }
  double fsec = std::strtod(start, &end);
  if (errno == ERANGE || end == start) {
    errno = 0;
    return 11;
  }

  t = ngpt::datetime<ngpt::microseconds>(
      ngpt::year(date[0]), ngpt::month(date[1]), ngpt::day_of_month(date[2]),
      ngpt::hours(date[3]), ngpt::minutes(date[4]),
      ngpt::microseconds(static_cast<long>(fsec * 1e6)));

  return 0;
}

/// Read in and resolve an Sp3c/d Velocity and ClockRate-of-Change Record. The
/// function expects that the next line to be read off from the input stream
/// is a Velocity and ClockRate-of-Change Record line.
///
/// @param[out] A 3character satellite id as recorded in the Sp3 file
/// @param[out] xv  X-component of satelite velocity, in dm/sec
/// @param[out] yv  y-component of satelite velocity, in dm/sec
/// @param[out] zv  Z-component of satelite velocity, in dm/sec
/// @param[out] cv  Clock rate-of-change in 10**-4 microseconds/second
/// @param[out] xstdv X-component std. deviation in 10**-4 mm/sec
/// @param[out] ystdv Y-component std. deviation in 10**-4 mm/sec
/// @param[out] zstdv Z-component std. deviation in 10**-4 mm/sec
/// @param[out] cstdv Clock std. deviation in 10**-4 psec/sec
/// @param[out] flag An Sp3Flag instance denoting the status of the resolved
///             fields. The flag is NOT reset (aka input flags will not be
///             touched). Any flags to be added, only affect position and clock
///             rate-of-change records (aka bad_abscent_velocity,
///             bad_abscent_clock_rate, has_vel_stddev, has_clk_rate_stdev).
/// @param[in] wsat If provided, then only resolve the data line if the given
///             SatelliteId wsat matches the one recorded in the line. If the
///             SatelliteId was not matched, the events bad_abscent_velocity
///             and bad_abscent_clock_rate are set.
/// @return Anything other than 0 denotes an error (note tha error codes must
///         be >0 and <10)
int ids::Sp3c::get_next_velocity(SatelliteId &sat, double &xv, double &yv,
                                 double &zv, double &cv, double &xstdv,
                                 double &ystdv, double &zstdv, double &cstdv,
                                 Sp3Flag &flag, const SatelliteId* wsat) noexcept {
  char line[MAX_RECORD_CHARS];
  char *end;
  const char *start;
  double dvec[4];

  __istream.getline(line, MAX_RECORD_CHARS);
  if (*line != 'V')
    return 1;

  std::memcpy(sat.id, line + 1, 3);

  if (wsat) {
    if (*wsat != sat) {
      flag.set(Sp3Event::bad_abscent_velocity | Sp3Event::bad_abscent_clock_rate);
      return 0;
    }
  }

  start = line + 4;
  for (int i = 0; i < 4; i++) {
    dvec[i] = std::strtod(start, &end);
    if (errno == ERANGE || end == start) {
      errno = 0;
      return 2 + i;
    }
    start += 14;
  }

  xv = dvec[0]; // dm/s
  yv = dvec[1];
  zv = dvec[2];
  cv = dvec[3]; // 10**-4 microseconds/second

  if (xv == 0e0 || (yv == 0e0 || zv == 0e0)) {
    flag.set(Sp3Event::bad_abscent_velocity);
  }

  if (cv >= SP3_MISSING_CLK_VALUE)
    flag.set(Sp3Event::bad_abscent_clock_rate);

  // std deviations (if any)
  int has_pos_stddev = false, has_clk_stddev = false;
  if (std::strlen(line) > 68) {
    if (*(line + 61) != ' ' || *(line + 62) != ' ') {
      int nn = std::strtol(line + 61, &end, 10);
      if (!nn || errno == ERANGE) {
        errno = 0;
        return 6;
      }
      xstdv = std::pow(fpb_pos__, nn); // 10**-4 mm/sec
      ++has_pos_stddev;
    }
    if (*(line + 64) != ' ' || *(line + 65) != ' ') {
      int nn = std::strtol(line + 64, &end, 10);
      if (!nn || errno == ERANGE) {
        errno = 0;
        return 6;
      }
      ystdv = std::pow(fpb_pos__, nn);
      ++has_pos_stddev;
    }
    if (*(line + 67) != ' ' || *(line + 68) != ' ') {
      int nn = std::strtol(line + 67, &end, 10);
      if (!nn || errno == ERANGE) {
        errno = 0;
        return 6;
      }
      zstdv = std::pow(fpb_pos__, nn);
      ++has_pos_stddev;
    }
  }
  if (has_pos_stddev == 3)
    flag.set(Sp3Event::has_vel_stddev);

  if (std::strlen(line) > 71) {
    if (*(line + 70) != ' ' || *(line + 71) != ' ' || *(line + 72) != ' ') {
      int nn = std::strtol(line + 70, &end, 10);
      if (!nn || errno == ERANGE) {
        errno = 0;
        return 6;
      }
      cstdv = std::pow(fpb_clk__, nn); // 10**-4 psec/sec
      ++has_clk_stddev;
    }
  }
  if (has_clk_stddev == 1)
    flag.set(Sp3Event::has_clk_rate_stdev);

  return 0;
}

/// Read in and resolve an Sp3c/d Position and Clock Record. The function
/// expects that the next line to be read off from the input stream is a
/// Position and Clock Record line.
///
/// @param[out] A 3character satellite id as recorded in the Sp3 file
/// @param[out] xkm X-component of satelite position, in km
/// @param[out] ykm y-component of satelite position, in km
/// @param[out] zkm Z-component of satelite position, in km
/// @param[out] clk Clock correction in microsec
/// @param[out] xstdv X-component std. deviation in mm
/// @param[out] ystdv Y-component std. deviation in mm
/// @param[out] zstdv Z-component std. deviation in mm
/// @param[out] cstdv Clock std. deviation in psec
/// @param[out] flag An Sp3Flag instance denoting the status of the resolved
///             fields. Note that the flag will be reset at the function call
/// @param[in] wsat If provided, then only resolve the data line if the given
///             SatelliteId wsat matches the one recorded in the line. If the
///             SatelliteId was not matched, the events bad_abscent_position
///             and bad_abscent_clock are set.
/// @return Anything other than 0 denotes an error (note tha error codes must
///         be >0 and <10)
int ids::Sp3c::get_next_position(SatelliteId &sat, double &xkm, double &ykm,
                                 double &zkm, double &clk, double &xstdv,
                                 double &ystdv, double &zstdv, double &cstdv,
                                 Sp3Flag &flag, const SatelliteId* wsat) noexcept {
  char line[MAX_RECORD_CHARS];
  char *end;
  const char *start;
  double dvec[4];
  
  flag.reset();

  __istream.getline(line, MAX_RECORD_CHARS);
  if (*line != 'P')
    return 1;

  std::memcpy(sat.id, line + 1, 3);
  
  if (wsat) {
    if (*wsat != sat) {
      flag.set(Sp3Event::bad_abscent_position | Sp3Event::bad_abscent_clock);
      return 0;
    }
  }

  start = line + 4;
  for (int i = 0; i < 4; i++) {
    dvec[i] = std::strtod(start, &end);
    if (errno == ERANGE || end == start) {
      errno = 0;
      return 2 + i;
    }
    start += 14;
  }

  xkm = dvec[0];
  ykm = dvec[1];
  zkm = dvec[2];
  clk = dvec[3];

  // std deviations (if any)
  int has_pos_stddev = false, has_clk_stddev = false;
  if (std::strlen(line) > 68) {
    if (*(line + 61) != ' ' || *(line + 62) != ' ') {
      int nn = std::strtol(line + 61, &end, 10);
      if (!nn || errno == ERANGE) {
        errno = 0;
        return 6;
      }
      xstdv = std::pow(fpb_pos__, nn);
      ++has_pos_stddev;
    }
    if (*(line + 64) != ' ' || *(line + 65) != ' ') {
      int nn = std::strtol(line + 64, &end, 10);
      if (!nn || errno == ERANGE) {
        errno = 0;
        return 6;
      }
      ystdv = std::pow(fpb_pos__, nn);
      ++has_pos_stddev;
    }
    if (*(line + 67) != ' ' || *(line + 68) != ' ') {
      int nn = std::strtol(line + 67, &end, 10);
      if (!nn || errno == ERANGE) {
        errno = 0;
        return 6;
      }
      zstdv = std::pow(fpb_pos__, nn);
      ++has_pos_stddev;
    }
  }
  if (has_pos_stddev == 3)
    flag.set(Sp3Event::has_pos_stddev);

  if (std::strlen(line) > 71) {
    if (*(line + 70) != ' ' || *(line + 71) != ' ' || *(line + 72) != ' ') {
      int nn = std::strtol(line + 70, &end, 10);
      if (!nn || errno == ERANGE) {
        errno = 0;
        return 6;
      }
      cstdv = std::pow(fpb_clk__, nn);
      ++has_clk_stddev;
    }
  }
  if (has_clk_stddev == 1)
    flag.set(Sp3Event::has_clk_stddev);

  if (xkm == 0e0 || (ykm == 0e0 || zkm == 0e0)) {
    flag.set(Sp3Event::bad_abscent_position);
  }

  if (clk >= SP3_MISSING_CLK_VALUE)
    flag.set(Sp3Event::bad_abscent_clock);
  if (line[74] == 'E')
    flag.set(Sp3Event::clock_event);
  if (line[75] == 'P')
    flag.set(Sp3Event::clock_prediction);
  if (line[78] == 'M')
    flag.set(Sp3Event::maneuver);
  if (line[79] == 'E')
    flag.set(Sp3Event::orbit_prediction);

  return 0;
}

/// @details Sp3c constructor, using a filename. The constructor will
///          initialize (set) the _filename attribute and also (try to)
///          open the input stream (i.e. _istream).
///          If the file is successefuly opened, the constructor will read
///          the header and assign info.
/// @param[in] filename  The filename of the Sp3 file
ids::Sp3c::Sp3c(const char *filename)
    : __filename(filename), __istream(filename, std::ios_base::in),
      /*__satsys(SATELLITE_SYSTEM::mixed),*/ __end_of_head(0) {
  int j;
  if ((j = read_header())) {
    if (__istream.is_open())
      __istream.close();
    throw std::runtime_error("[ERROR] Failed to read Sp3 header; Error Code: " +
                             std::to_string(j));
  }
}

/// @return -1: EOF encountered
///          0: All ok
///         >0: ERROR
int ids::Sp3c::get_next_data_block(SatelliteId satid, Sp3DataBlock& block) noexcept {
  char line[MAX_RECORD_CHARS];
  char c;
  int status;

  if (!__istream.good())
    return 1;

  /* possible following lines (three first chars):
   * 1. '*  ' i.e an epoch header
   * 2. 'PXX' i.e. a position & clock line, e.g. 'PG01 ....'
   * 3. 'EP ' i.e. position and clock correlation
   * 4. 'VXX' i.e. velocity line, e.g. 'VG01 ...'
   * 5. 'EV ' i.e. velocity correlation
   * 6. 'EOF' i.e. EOF
   */
  SatelliteId csatid;

  // following line should be an epoch header or 'EOF'
  c = __istream.peek();
  if (c == '*') {
    if ((status = resolve_epoch_line(block.t)))
      return status + 10;
  } else {
    __istream.getline(line, MAX_RECORD_CHARS);
    if (!std::strncmp(line, "EOF", 3)) {
      return -1;
    } else {
      return 100;
    }
  }

  // keep on reading reacords .....
  bool keep_reading = true;
  do {
    c = __istream.peek();
    if (c == '*') {
      keep_reading = false;
      break;
    } else if (c == 'P') {
      if ((status = get_next_position(satid, block.state[0], block.state[1], block.state[2],
                                      block.state[3], block.state_sdev[0], block.state_sdev[1],
                                      block.state_sdev[2], block.state_sdev[3], block.flag, &satid)))
        return status + 20;
    } else if (c == 'V') {
      if ((status = get_next_velocity(satid, block.state[4], block.state[5], block.state[6],
                                      block.state[7], block.state_sdev[4], block.state_sdev[5],
                                      block.state_sdev[6], block.state_sdev[7], block.flag, &satid)))
        return status + 30;
    } else {
      __istream.getline(line, MAX_RECORD_CHARS);
      if (!std::strncmp(line, "EOF", 3)) {
        keep_reading = false;
        return -1;
      } else if (!std::strncmp(line, "EP", 2)) {
        printf("[DEBUG] Ingoring Position Correlation Records ...\n");
      } else if (!std::strncmp(line, "EV", 2)) {
        printf("[DEBUG] Ingoring Velocity Correlation Records ...\n");
      } else {
        return 150;
      }
    }
  } while (keep_reading);

  return status;
}

#ifdef DEBUG
void ids::Sp3c::print_members() const noexcept {
  std::cout << "\nfilename     :" << __filename
            << "\nVersion      :" << version__
            << "\nStart Epoch  :" << ngpt::strftime_ymd_hms(start_epoch__)
            << "\n# Epochs     :" << num_epochs__
            << "\nCoordinate S :" << crd_sys__
            << "\nOrbit Type   :" << orb_type__
            << "\nAgency       :" << agency__
            << "\nTime System  :" << time_sys__
            << "\nInterval(sec):" << interval__.to_fractional_seconds();
}
#endif

ids::sp3_details::__Sp3FlagWrapper__ ids::operator|(ids::sp3_details::__Sp3FlagWrapper__ e1, ids::Sp3Event e2) noexcept {
  ids::sp3_details::__Sp3FlagWrapper__ wf;
  wf.bits_ = e1.bits_;
  wf.bits_ |= (1 << static_cast<ids::sp3_details::uitype>(e2));
  return wf;
}

ids::sp3_details::__Sp3FlagWrapper__ ids::operator|(ids::Sp3Event e1, ids::Sp3Event e2) noexcept {
  ids::sp3_details::__Sp3FlagWrapper__ wf;
  wf.bits_ |= (1 << static_cast<ids::sp3_details::uitype>(e1));
  wf.bits_ |= (1 << static_cast<ids::sp3_details::uitype>(e2));
  return wf;
}
