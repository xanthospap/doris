#include "doris_system_info.hpp"
#include <cassert>
#include <charconv>
#include <cstring>
#include <stdexcept>

namespace {
inline const char *skipws(const char *line) noexcept {
  const char *c = line;
  while (*c && *c == ' ')
    ++c;
  return c;
}
} // unnamed namespace

char dso::ObservationType_to_char(dso::ObservationType o) {
  switch (o) {
  case (ObservationType::phase):
    return 'L';
  case (ObservationType::pseudorange):
    return 'C';
  case (ObservationType::power_level):
    return 'W';
  case (ObservationType::frequency_offset):
    return 'F';
  case (ObservationType::ground_pressure):
    return 'P';
  case (ObservationType::ground_temperature):
    return 'T';
  case (ObservationType::ground_humidity):
    return 'H';
  default:
    throw std::runtime_error(
        "[ERROR] Cannot translate ObservationType to char");
  }
}

dso::ObservationType dso::char_to_observationType(char c) {
  switch (c) {
  case ('L'):
    return ObservationType::phase;
  case ('C'):
    return ObservationType::pseudorange;
  case ('W'):
    return ObservationType::power_level;
  case ('F'):
    return ObservationType::frequency_offset;
  case ('P'):
    return ObservationType::ground_pressure;
  case ('T'):
    return ObservationType::ground_temperature;
  case ('H'):
    return ObservationType::ground_humidity;
  default:
    throw std::runtime_error(
        "[ERROR] Cannot translate char to ObservationType");
  }
}

/// The function will check the type ObservationType instance, to see if
/// frequency information is needed for this ObservationType. That is, it will
/// return true if type is any of:
/// * ObservationType::phase, or
/// * ObservationType::pseudorange, or
/// * ObservationType::power_level
bool dso::observationType_has_frequency(ObservationType type) noexcept {
  switch (type) {
  case (ObservationType::phase):
  case (ObservationType::pseudorange):
  case (ObservationType::power_level):
    return true;
  case (ObservationType::frequency_offset):
  case (ObservationType::ground_pressure):
  case (ObservationType::ground_temperature):
  case (ObservationType::ground_humidity):
  default:
    return false;
  }
}

/// If the ObservationType is any of phase, pseudorange or power_level, then
/// the freq must be either 1 or 2. Otherwise (i.e. if type is not in phase,
/// pseudorange or power_level), the passed in frequency is ignored and set to
/// 0.
/// If the type is any of phase, pseudorange or power_level, and freq is not
/// in range [1,2], then the constructor will throw an std::runtime_error.
dso::ObservationCode::ObservationCode(ObservationType type, int_fast8_t freq)
    : m_type(type), m_freq(freq) {
  if (observationType_has_frequency(m_type)) {
    if (m_freq < 1 || m_freq > 2) {
      throw std::runtime_error("[ERROR] Invalid DORIS frequency number");
    }
  } else {
    m_freq = 0;
  }
}

/// The function will check the m_type member instance, to see if frequency
/// information is needed for this ObservationType. That is, it will return
/// true if m_type is any of:
/// * ObservationType::phase, or
/// * ObservationType::pseudorange, or
/// * ObservationType::power_level
bool dso::ObservationCode::has_frequency() const noexcept {
  return observationType_has_frequency(m_type);
}

int dso::BeaconStation::set_from_rinex_line(const char *line) noexcept {
  int sz;
  if ((sz = std::strlen(line)) < 60)
    return 1;
  std::memcpy(m_internal_code, line, sizeof m_internal_code);
  std::memcpy(m_station_id, line + 5, sizeof m_station_id);
  std::memcpy(m_station_name, line + 10, sizeof m_station_name);
  std::memcpy(m_station_domes, line + 40, sizeof m_station_domes);

  int error = 0;
  auto cerr = std::from_chars(skipws(line + 50), line + sz, m_type);
  if (cerr.ec != std::errc{})
    ++error;
  cerr = std::from_chars(skipws(line + 52), line + sz, m_type);
  if (cerr.ec != std::errc{})
    ++error;

  // remove trailing whitespace characters from stations name
  char *s = m_station_name + (sizeof m_station_name) - 1;
  while (*s == ' ' && (s - m_station_name) > 0)
    *s-- = '\0';
  s = m_station_domes + (sizeof m_station_domes) - 1;
  while (*s == ' ' && (s - m_station_domes) > 0)
    *s-- = '\0';

  return error;
}

char *dso::BeaconStation::to_str(char *buffer) const noexcept {
  std::sprintf(buffer, "[%.3s/%.4s_%.9s](%s)", m_internal_code, m_station_id,
               m_station_domes, m_station_name);
  return buffer;
}

char *dso::ObservationCode::to_str(char *buffer) const noexcept {
  assert(m_freq < 10);
  std::sprintf(buffer, "%c%1d", ObservationType_to_char(m_type), m_freq);
  return buffer;
}

dso::GroundAntennaType dso::BeaconStation::type() const {
  switch (m_station_id[3]) {
  case ('A'):
    return GroundAntennaType::Alcatel;
  case ('B'):
    return GroundAntennaType::Starec_B;
  case ('C'):
    return GroundAntennaType::Starec_C;
  default:
    throw std::runtime_error(
        "[ERROR] Cannot match beacon mnemonic to antenna type");
  }
}
