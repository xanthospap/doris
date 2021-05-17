#include "doris_rinex.hpp"
#include "ggdatetime/datetime_read.hpp"
#include <stdexcept>

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

int ids::DorisObsRinex::resolve_data_epoch(
    const char *line, ids::RinexDataRecordHeader &hdr) noexcept {
  if (*line != '>')
    return 1; // must start with '>' character

  char *end;
  int status = 0;

  try {
    hdr.m_epoch = ngpt::strptime_ymd_hms<ngpt::nanoseconds>(line + 2, &end);
  } catch (std::exception &e) {
    status = status ? status : 2;
  }

  hdr.m_flag = static_cast<int_fast8_t>(std::strtol(line + 33, &end, 10));
  if (end == line + 33 || errno) {
    errno = 0;
    status = status ? status : 3;
  }

  hdr.m_num_stations =
      static_cast<int_fast16_t>(std::strtol(line + 34, &end, 10));
  if (end == line + 34 || errno) {
    errno = 0;
    status = status ? status : 4;
  }

  bool has_clock_offset = false;
  for (int i = 0; i < 13; i++) {
    if (*(line + 42 + i) != ' ') {
      has_clock_offset = true;
      break;
    }
  }
  if (has_clock_offset) {
    hdr.m_clock_offset = std::strtod(line + 42, &end);
    if (errno || end == line + 42) {
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
