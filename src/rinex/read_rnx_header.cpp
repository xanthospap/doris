#include "datetime/datetime_read.hpp"
#include "doris_rinex.hpp"
#include <algorithm>
#include <cassert>
#include <cstring>

/// Max header lines.
/// constexpr int MAX_HEADER_LINES{1000};
namespace {
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
inline const char *skipws(const char *line) noexcept {
  const char *c = line;
  while (*c && *c == ' ')
    ++c;
  return c;
}
} // unnamed namespace

/// The instance's stream must be open and in good state. If it is not placed
/// at the top of the file, it will be rewinded to the top.
///
int dso::DorisObsRinex::read_header() noexcept {
  if (!m_stream.is_open() || !m_stream.good())
    return -1;
  if (m_stream.tellg())
    m_stream.seekg(0);

  char line[MAX_HEADER_CHARS];
  const char *start = line;
  int num_stations = 0;
  int num_ref_stations = 0;
  int obs_types_num = 0;
  int tmp_sz;

  // first line; RINEX VERSION / TYPE (get version, type and system)
  m_stream.getline(line, MAX_HEADER_CHARS);
  if (std::strncmp(line + 60, "RINEX VERSION / TYPE", 20))
    return 10;

  auto cres = std::from_chars(skipws(line), line + MAX_HEADER_CHARS, m_version);
  if (!(cres.ec == std::errc{}) || (line[20] != 'O' || line[40] != 'D'))
    return 11;

  // second line; PGM / RUN BY / DATE
  m_stream.getline(line, MAX_HEADER_CHARS);
  if (std::strncmp(line + 60, "PGM / RUN BY / DATE", 19))
    return 20;

  // read on untill EOH
  int error = 0;
  while (m_stream.getline(line, MAX_HEADER_CHARS) && !error) {
    if (!std::strncmp(line + 60, "SATELLITE NAME", 14)) {
      // SATELLITE NAME; get m_satellite_name (err. code 30)
      tmp_sz = count_length_reverse(line, 59);
      if (!tmp_sz || tmp_sz >= 59) {
        error = 30;
      } else {
        std::memcpy(m_satellite_name, line, tmp_sz);
        m_satellite_name[tmp_sz] = '\0';
      }

    } else if (!std::strncmp(line + 60, "COSPAR NUMBER", 13)) {
      // COSPAR NUMBER; get m_cospar_number (err. code 40)
      tmp_sz = count_length_reverse(line, 59);
      if (!tmp_sz || tmp_sz > 19) {
        error = 40;
      } else {
        std::memcpy(m_cospar_number, line, tmp_sz);
        m_cospar_number[tmp_sz] = '\0';
      }

    } else if (!std::strncmp(line + 60, "MARKER TYPE", 11)) {
      // MARKER TYPE; check that the field is "SPACEBORNE" (err. code 50)
      if (std::strncmp(line, "SPACEBORNE", 10))
        error = 50;

    } else if (!std::strncmp(line + 60, "OBSERVER / AGENCY", 17)) {
      // OBSERVER / AGENCY; currently ingored .... (err. code 60)
      ;

    } else if (!std::strncmp(line + 60, "REC # / TYPE / VERS", 19)) {
      // REC # / TYPE / VERS; get m_rec_chain, m_rec_type and m_rec_version
      // (err. code 70)
      tmp_sz = count_length_reverse(line, 19);
      if (tmp_sz < 1 || tmp_sz >= 19)
        error = 71;
      std::memcpy(m_rec_chain, line, tmp_sz);
      m_rec_chain[tmp_sz] = '\0';
      tmp_sz = count_length_reverse(line + 20, 19);
      if (tmp_sz < 1 || tmp_sz >= 19)
        error = 72;
      std::memcpy(m_rec_type, line + 20, tmp_sz);
      m_rec_type[tmp_sz] = '\0';
      tmp_sz = count_length_reverse(line + 40, 19);
      if (tmp_sz < 1 || tmp_sz >= 19)
        error = 73;
      std::memcpy(m_rec_version, line + 40, tmp_sz);
      m_rec_version[tmp_sz] = '\0';

    } else if (!std::strncmp(line + 60, "ANT # / TYPE", 12)) {
      // ANT # / TYPE; get and validate m_antenna_number and m_antenna_type
      // (err. code 80)
      tmp_sz = count_length_reverse(line, 19);
      if (tmp_sz < 1 || tmp_sz >= 19)
        error = 81;
      std::memcpy(m_antenna_number, line, tmp_sz);
      m_antenna_number[tmp_sz] = '\0';
      tmp_sz = count_length_reverse(line + 20, 19);
      if (tmp_sz < 1 || tmp_sz >= 19)
        error = 82;
      std::memcpy(m_antenna_type, line + 20, tmp_sz);
      m_antenna_type[tmp_sz] = '\0';
      if (std::strcmp(m_antenna_number, "DORIS") ||
          std::strcmp(m_antenna_type, "STAREC"))
        error = 84;

    } else if (!std::strncmp(line + 60, "APPROX POSITION XYZ", 19)) {
      // errno = 0;
      // APPROX POSITION XYZ; get m_approx_position (err. code 90)
      start = line;
      for (int i = 0; i < 3; i++) {
        cres = std::from_chars(skipws(start), line + 60, m_approx_position[i]);
        if (cres.ec != std::errc{}) {
          error = 90 + i;
        }
        start = cres.ptr;
      }

    } else if (!std::strncmp(line + 60, "CENTER OF MASS: XYZ", 19)) {
      // CENTER OF MASS: XYZ; get m_center_mass (err. code 100)
      start = line;
      for (int i = 0; i < 3; i++) {
        cres = std::from_chars(skipws(start), line + 60, m_center_mass[i]);
        if (cres.ec != std::errc{}) {
          error = 100 + i;
        }
        start = cres.ptr;
      }

    } else if (!std::strncmp(line + 60, "SYS / # / OBS TYPES", 19)) {
      // SYS / # / OBS TYPES; get/fill m_obs_codes (err. code 110)
      if (*line != 'D')
        return 111;
      cres = std::from_chars(skipws(line + 1), line + 60, obs_types_num);
      if (cres.ec != std::errc{}) {
        error = 112;
      }
      assert(obs_types_num < 13);
      if (!m_obs_codes.empty())
        m_obs_codes.clear();
      for (int code = 0; code < obs_types_num; code++) {
        char *s = line + 6 + code * 4;
        while (*s == ' ') // we need to skip leading whitespaces ...
          ++s;
        try {
          ObservationType type(char_to_observationType(*s)); // may throw!
          int freq = 0;
          if (observationType_has_frequency(type)) {
            cres = std::from_chars(skipws(s + 1), line + 60, freq);
            if (cres.ec != std::errc{}) {
              error = 113;
            }
          }
          m_obs_codes.emplace_back(type, freq); // may throw!
        } catch (std::exception &e) {
          error = 114;
        }
      }

    } else if (!std::strncmp(line + 60, "TIME OF FIRST OBS", 17)) {
      // TIME OF FIRST OBS; get m_time_of_first_obs and validate time system
      // (err. code 120)
      if (std::strncmp(line + 48, "DOR", 3))
        error = 121;
      try {
        m_time_of_first_obs = dso::strptime_ymd_hms<dso::nanoseconds>(line);
      } catch (std::exception &e) {
        error = 122;
      }

    } else if (!std::strncmp(line + 60, "SYS / DCBS APPLIED", 18)) {
      // SYS / DCBS APPLIED; not yet handled! (err. code 130)
      error = 130;

    } else if (!std::strncmp(line + 60, "SYS / SCALE FACTOR", 18)) {
      // SYS / SCALE FACTOR; get/fill m_obs_scale_factors
      if (*line != 'D')
        error = 141;
      if (!m_obs_scale_factors.empty())
        m_obs_scale_factors.clear();
      // create a vector of 1's with a size equal to m_obs_codes. the two
      // vectors will have a one-to-one correspondance
      m_obs_scale_factors = std::vector<int>(m_obs_codes.size(), 1);
      int factor;
      cres = std::from_chars(skipws(line + 2), line + 60, factor);
      if (cres.ec != std::errc{}) {
        error = 142;
      }
      int num_obs;
      if (line[8] == line[9] && line[8] == ' ') { // blank means all ....
        num_obs = m_obs_scale_factors.size();
      } else {
        // scale factors can be defined for any number of ObservationCodes, in
        // the range [0, m_obs_codes.size()]
        cres = std::from_chars(skipws(line + 8), line + 60, num_obs);
        if (cres.ec != std::errc{}) {
          error = 143;
        }
      }
      for (int code = 0; code < num_obs; code++) {
        char *s = line + 10 + code * 4;
        while (*s == ' ') // ignore leading whitespaces ...
          ++s;
        try {
          ObservationType type(char_to_observationType(*s));
          int freq = 0;
          if (observationType_has_frequency(type)) {
            cres = std::from_chars(skipws(s + 1), line + 60, freq);
            if (cres.ec != std::errc{}) {
              error = 144;
            }
          }
          ObservationCode tmp(type, freq);
          // get the index of the corresponding ObservationCode in the
          // m_obs_codes vector. It MUST be there ....
          auto it = std::find_if(
              m_obs_codes.cbegin(), m_obs_codes.cend(),
              [tmp](const ObservationCode &o) { return o == tmp; });
          if (it == m_obs_codes.cend())
            error = 145;
          m_obs_scale_factors[std::distance(m_obs_codes.cbegin(), it)] = factor;
        } catch (std::exception &e) {
          error = 146;
        }
      }

    } else if (!std::strncmp(line + 60, "L2 / L1 DATE OFFSET", 19)) {
      // L2 / L1 DATE OFFSET; get m_l12_date_offset (err. code 150)
      if (*line != 'D')
        error = 151;
      cres = std::from_chars(skipws(line + 3), line + 60, m_l12_date_offset);
      if (cres.ec != std::errc{}) {
        error = 152;
      }

    } else if (!std::strncmp(line + 60, "# OF STATIONS", 13)) {
      // # OF STATIONS; reserve size for m_stations (err. code 160)
      cres = std::from_chars(skipws(line), line + 60, num_stations);
      if (cres.ec != std::errc{}) {
        error = 161;
      }
      if ((int)m_stations.capacity() < num_stations)
        m_stations.reserve(num_stations);

    } else if (!std::strncmp(line + 60, "STATION REFERENCE", 17)) {
      // STATION REFERENCE; fill m_stations (err. code 170)
      m_stations.emplace_back(dso::BeaconStation{});
      if (m_stations[m_stations.size() - 1].set_from_rinex_line(line)) {
        error = 171;
      }

    } else if (!std::strncmp(line + 60, "# TIME REF STATIONS", 19)) {
      // # TIME REF STATIONS; reserve size for m_ref_stations (err. code 180)
      cres = std::from_chars(skipws(line), line + 60, num_ref_stations);
      if (cres.ec != std::errc{}) {
        error = 181;
      }
      if ((int)m_ref_stations.capacity() < num_ref_stations)
        m_ref_stations.reserve(num_ref_stations);

    } else if (!std::strncmp(line + 60, "TIME REF STATION", 16)) {
      // TIME REF STATION; collect time reference station. make sure reference
      // station is already in the m_stations vector (err. code 190)
      TimeReferenceStation refsta;
      std::memcpy(refsta.m_station_code, line, sizeof refsta.m_station_code);
      cres = std::from_chars(skipws(line + 5), line + 60, refsta.m_bias);
      if (cres.ec != std::errc{}) {
        error = 191;
      }
      start = cres.ptr;
      cres = std::from_chars(skipws(start), line + 60, refsta.m_shift);
      if (cres.ec != std::errc{}) {
        error = 192;
      }
      // find station/beacon in the m_stations vector
      if (auto it = std::find_if(m_stations.cbegin(), m_stations.cend(),
                                 [&refsta](const BeaconStation &s) {
                                   return !std::strncmp(
                                       s.m_internal_code, refsta.m_station_code,
                                       sizeof refsta.m_station_code);
                                 });
          it == m_stations.cend())
        error = 193;
      m_ref_stations.emplace_back(refsta);

    } else if (!std::strncmp(line + 60, "TIME REF STAT DATE",
                             18)) { // Code: 200
      // TIME REF STAT DATE; get m_time_ref_stat (err. code 200)
      try {
        m_time_ref_stat = dso::strptime_ymd_hms<dso::nanoseconds>(line);
      } catch (std::exception &e) {
        error = 201;
      }

    } else if (!std::strncmp(line + 60, "RCV CLOCK OFFS APPL",
                             19)) { // Code: 210

      int rcv_clock_offs_appl_int = -99;
      cres = std::from_chars(skipws(line), line + 60, rcv_clock_offs_appl_int);
      if ((cres.ec != std::errc{}) ||
          (rcv_clock_offs_appl_int < 0 || rcv_clock_offs_appl_int > 1)) {
        error = 210;
      }
      rcv_clock_offs_appl = rcv_clock_offs_appl_int;

    } else if (!std::strncmp(line + 60, "END OF HEADER", 13)) {
      // END OF HEADER; all done! break
      m_end_of_head = m_stream.tellg();
      break;

    } else {
      printf("[DEBUG] Ignoring header line [%s] (traceback: %s)\n", line,
             __func__);
    }

    if (error) {
      fprintf(
          stderr,
          "[ERROR] Error while reading RINEX header, code=%d (traceback: %s)\n",
          error, __func__);
    }
  }

  // check stream state
  if (!m_stream.good() || !m_end_of_head)
    return -2;

  if (error)
    return -3;

  // final checks on collected info
  if (m_obs_codes.size() != m_obs_scale_factors.size() ||
      obs_types_num != static_cast<int>(m_obs_codes.size()))
    return -4;
  if (num_stations != static_cast<int>(m_stations.size()) ||
      m_ref_stations.size() >= m_stations.size())
    return -4;

  return 0;
}
