#include "sp3.hpp"
#include <stdexcept>
#ifdef DEBUG
#include <iostream>
#include "ggdatetime/datetime_write.hpp"
#endif

/// Max record characters (for a navigation data block)
constexpr int MAX_RECORD_CHARS{128};

/// Bad or absent clock values areto be set to 999999.999999.  The six integer 
/// nines are required, whereasthe fractional part nines are optional.
constexpr double SP3_MISSING_CLK_VALUE {999999.e0};

bool substr_is_empty(const char* str, std::size_t count) noexcept {
  std::size_t idx = 0;
  while (idx<count && str[idx]==' ') ++idx;
  return idx == count;
}

/// @brief Resolve an Epoch Header Record line
/// @param[in] line An Epoch Header Record to be resolved
/// @param[out] t The epoch resolved from the input line
/// @return Anything other than 0 denotes an error
int resolve_epoch_line(const char* line, ngpt::datetime<ngpt::microseconds>& t) noexcept {
  if (line[0] != '*' || line[1] != ' ')
    return 2;
  
  int date[5];
  char *end, c;
  const char *start;
  
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

int ids::Sp3c::get_next_velocity(SatelliteId& sat, double& xv, double& yv, double& zv, double& cv, 
double& xstdv, double& ystdv, double& zstdv, double& cstdv, Sp3Flag& flag) noexcept {
  char line[MAX_RECORD_CHARS];
  char *end, c;
  const char* start;
  int status;
  double dvec[4];
  
  __istream.getline(line, MAX_RECORD_CHARS);
  if (*line != 'V')
    return 1;

  std::memcpy(sat.id, line+1, 3);
  
  start = line + 4;
  for (int i = 0; i < 4; i++) {
    dvec[i] = std::strtod(start, &end);
    if (errno == ERANGE || end == start) {
      errno = 0;
      return 5 + i;
    }
    start += 14;
  }

  xv = dvec[0]; // dm/s
  yv = dvec[1];
  zv = dvec[2];
  cv = dvec[3]; // 10**-4 microseconds/second
  
  if ( xv=0e0 || (yv=0e0 || zv==0e0) ) {
    flag.set(Sp3Event::bad_abscent_velocity);
  }

  if (cv >= SP3_MISSING_CLK_VALUE)
    flag.set(Sp3Event::bad_abscent_clock_rate);

  // std deviations (if any)
  int has_pos_stddev=false, has_clk_stddev=false;
  if (*(line+61)!=' ' || *(line+62)!=' ') {
    int nn = std::strtol(line + 61, &end, 10);
    if (!nn || errno == ERANGE) {
      errno = 0;
      return 15;
    }
    xstdv = std::pow(fpb_pos__, nn);
    ++has_pos_stddev;
  }
  if (*(line+64)!=' ' || *(line+65)!=' ') {
    int nn = std::strtol(line + 64, &end, 10);
    if (!nn || errno == ERANGE) {
      errno = 0;
      return 15;
    }
    ystdv = std::pow(fpb_pos__, nn);
    ++has_pos_stddev;
  }
  if (*(line+67)!=' ' || *(line+68)!=' ') {
    int nn = std::strtol(line + 67, &end, 10);
    if (!nn || errno == ERANGE) {
      errno = 0;
      return 15;
    }
    zstdv = std::pow(fpb_pos__, nn);
    ++has_pos_stddev;
  }
  if (has_pos_stddev==3) flag.set(Sp3Event::has_pos_stddev);


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
///             fields
/// @return Anything other than 0 denotes an error
int ids::Sp3c::get_next_position(SatelliteId& sat, double& xkm, double& ykm, double& zkm, double& clk, 
double& xstdv, double& ystdv, double& zstdv, double& cstdv, Sp3Flag& flag) noexcept {
  char line[MAX_RECORD_CHARS];
  char *end, c;
  const char* start;
  int status;
  double dvec[4];
  
  __istream.getline(line, MAX_RECORD_CHARS);
  if (*line != 'P')
    return 1;

  std::memcpy(sat.id, line+1, 3);

  start = line + 4;
  for (int i = 0; i < 4; i++) {
    dvec[i] = std::strtod(start, &end);
    if (errno == ERANGE || end == start) {
      errno = 0;
      return 5 + i;
    }
    start += 14;
  }

  xkm = dvec[0];
  ykm = dvec[1];
  zkm = dvec[2];
  clk = dvec[3];

  flag.reset();
  
  // std deviations (if any)
  int has_pos_stddev=false, has_clk_stddev=false;
  if (*(line+61)!=' ' || *(line+62)!=' ') {
    int nn = std::strtol(line + 61, &end, 10);
    if (!nn || errno == ERANGE) {
      errno = 0;
      return 15;
    }
    xstdv = std::pow(fpb_pos__, nn);
    ++has_pos_stddev;
  }
  if (*(line+64)!=' ' || *(line+65)!=' ') {
    int nn = std::strtol(line + 64, &end, 10);
    if (!nn || errno == ERANGE) {
      errno = 0;
      return 15;
    }
    ystdv = std::pow(fpb_pos__, nn);
    ++has_pos_stddev;
  }
  if (*(line+67)!=' ' || *(line+68)!=' ') {
    int nn = std::strtol(line + 67, &end, 10);
    if (!nn || errno == ERANGE) {
      errno = 0;
      return 15;
    }
    zstdv = std::pow(fpb_pos__, nn);
    ++has_pos_stddev;
  }
  if (has_pos_stddev==3) flag.set(Sp3Event::has_pos_stddev);

  if (*(line+70)!=' ' || *(line+71)!=' ' || *(line+72)!=' ') {
    int nn = std::strtol(line + 70, &end, 10);
    if (!nn || errno == ERANGE) {
      errno = 0;
      return 15;
    }
    cstdv = std::pow(fpb_clk__, nn);
    ++has_clk_stddev;
  }
  if (has_clk_stddev==1) flag.set(Sp3Event::has_clk_stddev);
  
  if ( xkm=0e0 || (ykm=0e0 || zkm==0e0) ) {
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

struct Sp3DataBlock {
  ngpt::datetime<ngpt::microseconds> t__;
};

/// @return -1: EOF encountered
///          0: All ok
//.         >0: ERROR
int ids::Sp3c::get_next_data_block() noexcept {
  char line[MAX_RECORD_CHARS];
  char *end, *start, c;
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
  
  // resolve epoch header line
  __istream.getline(line, MAX_RECORD_CHARS);
  if (!std::strncmp(line, "EOF", 3)) { /* EOF */
    return -1;
  }
  if ((status = resolve_epoch_line(line, t))) {
    return status;
  }

  // peek next line's char to see what kind of line we are reading
  do {
    c = __istream.peek();
    if (c != '*') {
      __istream.getline(line, MAX_RECORD_CHARS);
      if (*line == 'P') { /* position and clock line */
        if ((j = get_next_position()))
          return j + 20;
      } else { /* any other line is irrelevant */
        ;
      }
    } else { /* next line is header line */
      keep_reading = false;
    }
  } while (keep_reading);

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
