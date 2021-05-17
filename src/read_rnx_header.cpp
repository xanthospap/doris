#include "doris_rinex.hpp"
#include "ggdatetime/datetime_read.hpp"
#include <algorithm>
#include <cassert>
#include <cstring>

/// No header line can have more than 80 chars.
constexpr int MAX_HEADER_CHARS{81};

/// Max header lines.
/// constexpr int MAX_HEADER_LINES{1000};

/// Starting from position from in the given string str, find the first non-
/// whitespace character, and return the length of str from start to this
/// position.
/// Examples:
/// str="foobar                 ", from=22, return value=6
/// str="foobar  barfoo         ", from=22, return value=14
/// str="   bar  barfoo         ", from=22, return value=14
/// str="                       ", from=22, return value=0
/// str=" a                     ", from=22, return value=2
int count_length_reverse(const char *str, int from) noexcept {
  while (*(str + from) == ' ' && from >= 0)
    --from;
  return ++from;
}

/// The instance's stream must be open and in good state. If it is not placed
/// at the top of the file, it will be rewinded to the top.
///
int ids::DorisObsRinex::read_header() noexcept {
  if (!m_stream.is_open() || !m_stream.good())
    return -1;
  if (m_stream.tellg())
    m_stream.seekg(0);

  char line[MAX_HEADER_CHARS];
  char *str_end, *end;
  int num_stations = 0;
  int num_ref_stations = 0;
  int obs_types_num = 0;
  int tmp_sz;

  // first line; RINEX VERSION / TYPE (get version, type and system)
  m_stream.getline(line, MAX_HEADER_CHARS);
  if (std::strncmp(line + 60, "RINEX VERSION / TYPE", 20))
    return 10;
  m_version = std::strtof(line, &str_end);
  if (str_end == line || (line[20] != 'O' || line[40] != 'D'))
    return 11;

  // second line; PGM / RUN BY / DATE
  m_stream.getline(line, MAX_HEADER_CHARS);
  if (std::strncmp(line + 60, "PGM / RUN BY / DATE", 19))
    return 20;

  // read on untill EOH
  while (m_stream && m_stream.getline(line, MAX_HEADER_CHARS)) {
    if (!std::strncmp(line + 60, "SATELLITE NAME", 14)) {
      // SATELLITE NAME; get m_satellite_name (err. code 30)
      tmp_sz = count_length_reverse(line, 59);
      if (!tmp_sz || tmp_sz >= 59)
        return 30;
      std::memcpy(m_satellite_name, line, tmp_sz);
      m_satellite_name[tmp_sz] = '\0';
    } else if (!std::strncmp(line + 60, "COSPAR NUMBER", 13)) {
      // COSPAR NUMBER; get m_cospar_number (err. code 40)
      tmp_sz = count_length_reverse(line, 59);
      if (!tmp_sz || tmp_sz > 19)
        return 40;
      std::memcpy(m_cospar_number, line, tmp_sz);
      m_cospar_number[tmp_sz] = '\0';
    } else if (!std::strncmp(line + 60, "MARKER TYPE", 11)) {
      // MARKER TYPE; check that the field is "SPACEBORNE" (err. code 50)
      if (std::strncmp(line, "SPACEBORNE", 10))
        return 50;
    } else if (!std::strncmp(line + 60, "OBSERVER / AGENCY", 17)) {
      // OBSERVER / AGENCY; currently ingored .... (err. code 60)
      ;
    } else if (!std::strncmp(line + 60, "REC # / TYPE / VERS", 19)) {
      // REC # / TYPE / VERS; get m_rec_chain, m_rec_type and m_rec_version
      // (err. code 70)
      tmp_sz = count_length_reverse(line, 19);
      if (tmp_sz < 1 || tmp_sz >= 19)
        return 71;
      std::memcpy(m_rec_chain, line, tmp_sz);
      m_rec_chain[tmp_sz] = '\0';
      tmp_sz = count_length_reverse(line + 20, 19);
      if (tmp_sz < 1 || tmp_sz >= 19)
        return 72;
      std::memcpy(m_rec_type, line + 20, tmp_sz);
      m_rec_type[tmp_sz] = '\0';
      tmp_sz = count_length_reverse(line + 40, 19);
      if (tmp_sz < 1 || tmp_sz >= 19)
        return 73;
      std::memcpy(m_rec_version, line + 40, tmp_sz);
      m_rec_version[tmp_sz] = '\0';
    } else if (!std::strncmp(line + 60, "ANT # / TYPE", 12)) {
      // ANT # / TYPE; get and validate m_antenna_number and m_antenna_type
      // (err. code 80)
      tmp_sz = count_length_reverse(line, 19);
      if (tmp_sz < 1 || tmp_sz >= 19)
        return 81;
      std::memcpy(m_antenna_number, line, tmp_sz);
      m_antenna_number[tmp_sz] = '\0';
      tmp_sz = count_length_reverse(line + 20, 19);
      if (tmp_sz < 1 || tmp_sz >= 19)
        return 82;
      std::memcpy(m_antenna_type, line + 20, tmp_sz);
      m_antenna_type[tmp_sz] = '\0';
      if (std::strcmp(m_antenna_number, "DORIS") ||
          std::strcmp(m_antenna_type, "STAREC"))
        return 84;
    } else if (!std::strncmp(line + 60, "APPROX POSITION XYZ", 19)) {
      // APPROX POSITION XYZ; get m_approx_position (err. code 90)
      char *start = line;
      for (int i = 0; i < 3; i++) {
        m_approx_position[i] = std::strtof(start, &end);
        if (start == end || errno) {
          errno = 0;
          return 90 + i;
        }
        start += 14;
      }
    } else if (!std::strncmp(line + 60, "CENTER OF MASS: XYZ", 19)) {
      // CENTER OF MASS: XYZ; get m_center_mass (err. code 100)
      char *start = line;
      for (int i = 0; i < 3; i++) {
        m_center_mass[i] = std::strtof(start, &end);
        if (start == end || errno) {
          errno = 0;
          return 100 + i;
        }
        start += 14;
      }
    } else if (!std::strncmp(line + 60, "SYS / # / OBS TYPES", 19)) {
      // SYS / # / OBS TYPES; get/fill m_obs_codes (err. code 110)
      if (*line != 'D')
        return 111;
      obs_types_num = std::strtol(line + 1, &end, 10);
      if ((errno || end == line + 1) || !obs_types_num) {
        errno = 0;
        return 112;
      }
      assert(obs_types_num < 13);
      if (!m_obs_codes.empty())
        m_obs_codes.clear();
      m_obs_codes.reserve(obs_types_num);
      for (int code = 0; code < obs_types_num; code++) {
        char *s = line + 6 + code * 4;
        while (*s == ' ') // we need to skip leading whitespaces ...
          ++s;
        try {
          ObservationType type(char_to_observationType(*s)); // may throw!
          int freq = 0;
          if (observationType_has_frequency(type)) {
            freq = std::strtol(s + 1, &end, 10);
            if (errno || end == s + 1) {
              errno = 0;
              return 113;
            }
          }
          m_obs_codes.emplace_back(type, freq); // may throw!
        } catch (std::exception &e) {
          return 114;
        }
      }
    } else if (!std::strncmp(line + 60, "TIME OF FIRST OBS", 17)) {
      // TIME OF FIRST OBS; get m_time_of_first_obs and validate time system
      // (err. code 120)
      if (std::strncmp(line + 48, "DOR", 3))
        return 121;
      try {
        m_time_of_first_obs = ngpt::strptime_ymd_hms<ngpt::nanoseconds>(line);
      } catch (std::exception &e) {
        return 122;
      }
    } else if (!std::strncmp(line + 60, "SYS / DCBS APPLIED", 18)) {
      // SYS / DCBS APPLIED; not yet handled! (err. code 130)
      return 130;
    } else if (!std::strncmp(line + 60, "SYS / SCALE FACTOR", 18)) {
      // SYS / SCALE FACTOR; get/fill m_obs_scale_factors
      if (*line != 'D')
        return 141;
      if (!m_obs_scale_factors.empty())
        m_obs_scale_factors.clear();
      // create a vector of 1's with a size equal to m_obs_codes. the two
      // vectors will have a one-to-one correspondance
      m_obs_scale_factors = std::vector<int>(m_obs_codes.size(), 1);
      int factor = std::strtol(line + 2, &end, 10);
      if (!factor || (errno || end == line + 2)) {
        errno = 0;
        return 142;
      }
      int num_obs;
      if (line[8] == line[9] && line[8] == ' ') { // blank means all ....
        num_obs = m_obs_scale_factors.size();
      } else {
        // scale factors can be defined for any number of ObservationCodes, in
        // the range [0, m_obs_codes.size()]
        num_obs = std::strtol(line + 8, &end, 10);
        if ((!num_obs || num_obs > (int)m_obs_scale_factors.size()) ||
            (errno || end == line + 8))
          return 143;
      }
      for (int code = 0; code < num_obs; code++) {
        char *s = line + 10 + code * 4;
        while (*s == ' ') // ignore leading whitespaces ...
          ++s;
        try {
          ObservationType type(char_to_observationType(*s));
          int freq = 0;
          if (observationType_has_frequency(type)) {
            freq = std::strtol(s + 1, &end, 10);
            if (!freq || (errno || end == s + 1)) {
              errno = 0;
              return 144;
            }
          }
          ObservationCode tmp(type, freq);
          // get the index of the corresponding ObservationCode in the
          // m_obs_codes vector. It MUST be there ....
          auto it = std::find_if(
              m_obs_codes.cbegin(), m_obs_codes.cend(),
              [tmp](const ObservationCode &o) { return o == tmp; });
          if (it == m_obs_codes.cend())
            return 145;
          m_obs_scale_factors[std::distance(m_obs_codes.cbegin(), it)] = factor;
        } catch (std::exception &e) {
          return 146;
        }
      }
    } else if (!std::strncmp(line + 60, "L2 / L1 DATE OFFSET", 19)) {
      // L2 / L1 DATE OFFSET; get m_l12_date_offset (err. code 150)
      if (*line != 'D')
        return 151;
      m_l12_date_offset = std::strtod(line + 3, &end);
      if (errno || end == line + 3) {
        errno = 0;
        return 152;
      }
    } else if (!std::strncmp(line + 60, "# OF STATIONS", 13)) {
      // # OF STATIONS; reserve size for m_stations (err. code 160)
      num_stations = std::strtol(line, &end, 10);
      if (!num_stations || (errno || end == line))
        return 161;
      m_stations.reserve(num_stations);
    } else if (!std::strncmp(line + 60, "STATION REFERENCE", 17)) {
      // STATION REFERENCE; fill m_stations (err. code 170)
      m_stations.emplace_back();
      if (m_stations[m_stations.size() - 1].set_from_rinex_line(line)) {
        return 171;
      }
    } else if (!std::strncmp(line + 60, "# TIME REF STATIONS", 19)) {
      // # TIME REF STATIONS; reserve size for m_ref_stations (err. code 180)
      num_ref_stations = std::strtol(line, &end, 10);
      if (!num_ref_stations || (errno || end == line)) {
        errno = 0;
        return 181;
      }
      m_ref_stations.reserve(num_ref_stations);
    } else if (!std::strncmp(line + 60, "TIME REF STATION", 16)) {
      // TIME REF STATION; collect time reference station. make sure reference
      // station is already in the m_stations vector (err. code 190)
      TimeReferenceStation refsta;
      std::memcpy(refsta.m_station_code, line, sizeof refsta.m_station_code);
      refsta.m_bias = std::strtod(line + 5, &end);
      if (errno || end == line + 5) {
        errno = 0;
        return 191;
      }
      refsta.m_shift = std::strtod(line + 5 + 14, &end);
      if (errno || end == line + 5 + 14) {
        errno = 0;
        return 192;
      }
      // find station/beacon in the m_stations vector
      if (auto it = std::find_if(m_stations.cbegin(), m_stations.cend(),
                                 [&refsta](const BeaconStation &s) {
                                   return !std::strncmp(
                                       s.m_internal_code, refsta.m_station_code,
                                       sizeof refsta.m_station_code);
                                 });
          it == m_stations.cend())
        return 193;
      m_ref_stations.emplace_back(refsta);
    } else if (!std::strncmp(line + 60, "TIME REF STAT DATE",
                             18)) { // Code: 200
      // TIME REF STAT DATE; get m_time_ref_stat (err. code 200)
      try {
        m_time_ref_stat = ngpt::strptime_ymd_hms<ngpt::nanoseconds>(line);
      } catch (std::exception &e) {
        return 201;
      }
    } else if (!std::strncmp(line + 60, "END OF HEADER", 13)) {
      // END OF HEADER; all done! break
      m_end_of_head = m_stream.tellg();
      break;
    } else {
      std::cout << "\n[DEBUG] Skipping RINEX header line:\n" << line;
    }
  }

  // check stream state
  if (!m_stream.good() || !m_end_of_head)
    return -2;

  // final checks on collected info
  if (m_obs_codes.size() != m_obs_scale_factors.size() ||
      obs_types_num != static_cast<int>(m_obs_codes.size()))
    return -3;
  if (num_stations != static_cast<int>(m_stations.size()) ||
      m_ref_stations.size() >= m_stations.size())
    return -3;

  return 0;
}
