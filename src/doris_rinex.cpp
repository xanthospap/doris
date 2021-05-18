#include "doris_rinex.hpp"
#include "ggdatetime/datetime_read.hpp"
#include <stdexcept>

/// The constructor will try to:
/// 1. open the input file
/// 2. parse the header
/// If any of the above fails, then an std::runtime_error will be thrown.
ids::DorisObsRinex::DorisObsRinex(const char *fn)
    : m_filename(fn), m_stream(fn, std::ios_base::in) {
  int status = read_header();
  if (status) {
    std::cerr << "\n[ERROR] Failed reading RINEX header; filename : \"" << fn
              << "\", Error Code: " << status;
    throw std::runtime_error("[ERROR] Cannot read RINEX header");
  }
}

void ids::DorisObsRinex::skip_next_epoch(int num_stations,
                                         int lines_per_station) noexcept {
  constexpr int MAX_RECORD_CHARS = 124;
  char line[MAX_RECORD_CHARS];
  for (int i = 0; i < num_stations * lines_per_station; i++) {
    m_stream.getline(line, MAX_RECORD_CHARS);
  }
  return;
}

void ids::DorisObsRinex::read() {
  m_stream.seekg(m_end_of_head);
  constexpr int MAX_RECORD_CHARS = 124;
  char line[MAX_RECORD_CHARS];
  RinexDataRecordHeader hdr;

  while (m_stream && m_stream.getline(line, MAX_RECORD_CHARS)) {
    if (m_stream.eof())
      break;
    if (int status = resolve_data_epoch(line, hdr); status) {
      std::cerr << "\nFailed to read record header";
      std::cerr << "\nline: " << line;
      std::cerr << "\nError: " << status;
      return;
    }
    std::cout << "\nResolved line: \"" << line << "\"";
    std::cout << "\n\t" << hdr.m_clock_offset << ", " << hdr.m_num_stations
              << ", " << (int)hdr.m_flag << ", " << (int)hdr.m_clock_flag;
    skip_next_epoch(hdr.m_num_stations, this->lines_per_beacon());
  }

  return;
}

/// Example line:
/*
> 2020 01 01 01 41 53.279947800  0  4       -4.432841287 0 
  +-------------------------+-----------+-----------+----+
  | Record identifier : '>' |  A1       | start:  0 | 1
  +-------------------------+-----------+-----------+----+
  | Epoch :                 |           |
  | - year (4 digits)       | 1X,I4     | start:  1 | 5
  | - month,day,hour,min    | 4(1X,I2.2)| start:  6 | 12
  | - sec                   | F13.9     | start: 18 | 13
  +-------------------------+-----------+-----------+----+
  |Epoch flag               | 2X,I1     | start: 31 | 3
  |   0: OK                             |           |
  |   1: power failure between          |           |
  |      previous and current epoch     |           |
  |  >1: Special event                  |           |
  +-------------------------+-----------+-----------+----+
  |Number of stations       | I3        | start: 34 | 3
  |observed in current epoch|           |           |
  +-------------------------+-----------+-----------+----+
  |(reserved)               | 6X        | start: 37 | 6
  +-------------------------+-----------+-----------+----+
  | Receiver clock offset   | F13.9     | start: 43 | 13
  | (seconds, optional)     |           |           |
  +-------------------------+-----------+-----------+----+
  | Receiver clock offset   |           | start: 56 | 3
  | flag,                   | 1X,I1,1X  |           |
  |  - 1 if extrapolated,   |           |           |
  |  - 0 otherwise          |           | Max length of line = 59 chars
  +-------------------------+-----------+------------------------------
*/
int ids::DorisObsRinex::resolve_data_epoch(
    const char *line, ids::RinexDataRecordHeader &hdr) const noexcept {
  if (*line != '>')
    return 1; // must start with '>' character

  char *end;
  int status = 0;

  try {
    hdr.m_epoch = ngpt::strptime_ymd_hms<ngpt::nanoseconds>(line + 2, &end);
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
  hdr.m_num_stations =
         static_cast<int_fast16_t>(std::strtol(tbuf, &end, 10));
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
