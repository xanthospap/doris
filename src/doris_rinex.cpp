#include "doris_rinex.hpp"
#include "ggdatetime/datetime_read.hpp"
#include <cctype>
#include <stdexcept>
#ifdef DEBUG
#include "ggdatetime/datetime_write.hpp"
#endif

/// The constructor will try to:
/// 1. open the input file
/// 2. parse the header
/// If any of the above fails, then an std::runtime_error will be thrown.
ids::DorisObsRinex::DorisObsRinex(const char *fn)
    : m_filename(fn), m_stream(fn, std::ios_base::in) {
  int status = read_header();
  if (status) {
    fprintf(stderr, "[ERROR] Failed reading RINEX header for %s (error=%d) (traceback: %s)\n", fn, status, __func__);
    throw std::runtime_error("[ERROR] Cannot read RINEX header");
  }
  m_lines_per_beacon = lines_per_beacon();
}

void ids::DorisObsRinex::skip_data_block(
    const ids::RinexDataRecordHeader &hdr) noexcept {
  char line[MAX_RECORD_CHARS];
  for (int i = 0; i < hdr.m_num_stations * m_lines_per_beacon; i++) {
    m_stream.getline(line, MAX_RECORD_CHARS);
  }
  return;
}

#ifdef DEBUG
void ids::DorisObsRinex::read() {
  m_stream.seekg(m_end_of_head);
  char line[MAX_RECORD_CHARS];
  RinexDataRecordHeader hdr;
  std::vector<BeaconObservations> obsvec;
  dso::datetime<dso::nanoseconds> last_epoch;

  while (m_stream && m_stream.getline(line, MAX_RECORD_CHARS)) {
    if (m_stream.eof())
      break;
    if (int status = resolve_data_epoch(line, hdr); status) {
      std::cerr << "\nFailed to read record header";
      std::cerr << "\nline: " << line;
      std::cerr << "\nError: " << status;
      return;
    }
    hdr.apply_clock_offset();
    std::cout << "Read Epoch : "<<dso::strftime_ymd_hms(hdr.m_epoch) ;
    std::cout << " diff in seconds: " << hdr.m_epoch.delta_sec(last_epoch).as_underlying_type()/1000000000L << "\n";
    last_epoch = hdr.m_epoch;
    // std::cout << "\nResolved line: \"" << line << "\"";
    // std::cout << "\n\t" << hdr.m_clock_offset << ", " << hdr.m_num_stations
    //           << ", " << (int)hdr.m_flag << ", " << (int)hdr.m_clock_flag;
    // skip_next_epoch(hdr.m_num_stations, m_lines_per_beacon);
    if (int status = read_data_block(hdr, obsvec); status) {
      std::cerr << "\nFailed to read data record! error code " << status;
      return;
    }
  }

  return;
}
#endif

/// Example:
/// D21  -1339276.297 0   -263905.449 0 115619271.79511 115619357.25411      -119.750 7
///          -108.200 7      2149.340         530.000 1       -10.000 1        60.000 1
/// +---------------------------------------------------------------------------+
/// 
/// * A1,I2 aka [0-3)   Station number e.g. 'D02' Only at first line, else 3X
/// * F14.3 aka [3-17)  Observation
/// * I1    aka [17-18) flag-m1
/// * I2    aka [18-19) flag m2
/// if more than 5 measurements then repeat (instead of Station number we have 3X)
/// Missing observations are written as 0.0 or blanks.
/// maximum number of chars per (record) line: 3 + 5*16 = 83
/// @note Note that the resulting observable values are always scaled using the 
///       'SYS / SCALE FACTOR' header information (stored in m_obs_scale_factors
int ids::DorisObsRinex::read_data_block(
    ids::RinexDataRecordHeader &hdr,
    std::vector<BeaconObservations> &obsvec) noexcept {
  static char line[MAX_RECORD_CHARS];
  if (!obsvec.empty())
    obsvec.clear();

  // temporary buffer to hold fields to be resolved
  char buf[16] = {'\0'};
  char *end;
  double val;

  // loop for every beacon/station in the RinexDataRecordHeader ...
  for (int beacon = 0; beacon < hdr.m_num_stations; beacon++) {
    // append a new BeaconObservations instance for the current beacon ...
    obsvec.emplace_back(m_obs_codes.size());
    // and get an iterator to it (so that we set its values in-place)
    auto obsvec_it = obsvec.end() - 1;
    int curobs = 0;
    int curline = 0;
    // for every observation code descrbed in the RINEX header ...
    while (curobs < (int)m_obs_codes.size()) {
      // should we change/get the next line ?
      if (!(curobs % 5)) {
        m_stream.getline(line, MAX_RECORD_CHARS);
        // if this is the first data line for the beacon, get its code
        if (!curline) {
          if ((*line) != 'D')
            return 1;
          std::memcpy(obsvec_it->m_beacon_id, line, 3);
        }
        ++curline;
      }
      // collect measurements ...
      std::memcpy(buf, line + 3 + (curobs % 5) * 16, 14);
      // check if value is ommited (an ommited value is either left blank, or
      // is recorded as 0.0)
      bool buf_is_empty = true;
      end = buf;
      while (*end) {
        if (!isspace(*end++)) {
          buf_is_empty = false;
          break;
        }
      }
      char flagm1 = line[3 + (curobs % 5) * 16 + 14];
      char flagm2 = line[3 + (curobs % 5) * 16 + 15];
      val = (buf_is_empty) ? OBSERVATION_VALUE_MISSING : std::strtod(buf, &end);
      if (val == 0e0 || end == buf) {
        val = OBSERVATION_VALUE_MISSING;
        if (end == buf)
          return 2;
      }
      // push value to the current BeaconObservations instance (in-place)
      // WAIT! check if we have a scale factor for the observable (note that
      // the m_obs_scale_factors are in one-to-one correspondance with the
      // m_obs_codes vector. Hence, we can find the scale factor simply by the
      // index of the observable. If no scale factor for the observable exists,
      // then the m_obs_scale_factors should have an '1' in the corresponding
      // index
      val /= m_obs_scale_factors[curobs];
      obsvec_it->m_values.emplace_back(val, flagm1, flagm2);
      ++curobs;
    }
  }

  /*
  #ifdef DEBUG
    for (const auto& bcn: obsvec) {
      std::cout<<"\n\t->"<<bcn.to_string();
    }
  #endif
  */
  return 0;
}

/// Example line:
///
/// > 2020 01 01 01 41 53.279947800  0  4       -4.432841287 0
///   +-------------------------+-----------+-----------+----+
///   | Record identifier : '>' |  A1       | start:  0 | 1
///   +-------------------------+-----------+-----------+----+
///   | Epoch :                 |           |
///   | - year (4 digits)       | 1X,I4     | start:  1 | 5
///   | - month,day,hour,min    | 4(1X,I2.2)| start:  6 | 12
///   | - sec                   | F13.9     | start: 18 | 13
///   +-------------------------+-----------+-----------+----+
///   |Epoch flag               | 2X,I1     | start: 31 | 3
///   |   0: OK                             |           |
///   |   1: power failure between          |           |
///   |      previous and current epoch     |           |
///   |  >1: Special event                  |           |
///   +-------------------------+-----------+-----------+----+
///   |Number of stations       | I3        | start: 34 | 3
///   |observed in current epoch|           |           |
///   +-------------------------+-----------+-----------+----+
///   |(reserved)               | 6X        | start: 37 | 6
///   +-------------------------+-----------+-----------+----+
///   | Receiver clock offset   | F13.9     | start: 43 | 13
///   | (seconds, optional)     |           |           |
///   +-------------------------+-----------+-----------+----+
///   | Receiver clock offset   |           | start: 56 | 3
///   | flag,                   | 1X,I1,1X  |           |
///   |  - 1 if extrapolated,   |           |           |
///   |  - 0 otherwise          |           | Max length of line = 59 chars
///   +-------------------------+-----------+------------------------------
///
int ids::DorisObsRinex::resolve_data_epoch(
    const char *line, ids::RinexDataRecordHeader &hdr) const noexcept {
  if (*line != '>')
    return 1; // must start with '>' character

  char *end;
  int status = 0;

  try {
    hdr.m_epoch = dso::strptime_ymd_hms<dso::nanoseconds>(line + 2, &end);
  } catch (std::exception &e) {
    status = status ? status : 2;
  }

  // probably i am exagerating a bit here, but nevertheless ...
  // it could happen that the 'Epoch flag' and the 'Number of stations' fields
  // are joined in one big int (if number of stations is >=100). Hence, just to
  // be safe, we are moving the firlds in a temporary buffer and parse from
  // there
  char tbuf[4];
  std::memset(tbuf, '\0', sizeof tbuf);
  std::memcpy(tbuf, line + 31, 3);
  // hdr.m_flag = static_cast<int_fast8_t>(std::strtol(line + 31, &end, 10));
  hdr.m_flag = static_cast<int_fast8_t>(std::strtol(tbuf, &end, 10));
  if (end == tbuf || errno) {
    errno = 0;
    status = status ? status : 3;
  }

  std::memcpy(tbuf, line + 34, 3);
  // hdr.m_num_stations =
  //     static_cast<int_fast16_t>(std::strtol(line + 34, &end, 10));
  hdr.m_num_stations = static_cast<int_fast16_t>(std::strtol(tbuf, &end, 10));
  if (end == tbuf || errno) {
    errno = 0;
    status = status ? status : 4;
  }

  bool has_clock_offset = false;
  for (int i = 0; i < 13; i++) {
    if (*(line + 43 + i) != ' ') {
      has_clock_offset = true;
      break;
    }
  }
  if (has_clock_offset) {
    hdr.m_clock_offset = std::strtod(line + 43, &end);
    if (errno || end == line + 43) {
      errno = 0;
      status = status ? status : 5;
    }
  } else {
    hdr.m_clock_offset = RECEIVER_CLOCK_OFFSET_MISSING;
  }

  hdr.m_clock_flag = static_cast<int_fast8_t>(std::strtol(line + 56, &end, 10));
  if (errno || end == line + 56) {
    errno = 0;
    status = status ? status : 6;
  }

  return status;
}
