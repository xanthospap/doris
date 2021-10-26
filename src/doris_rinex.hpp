#ifndef __IDS_DORIS_RINEX_HPP__
#define __IDS_DORIS_RINEX_HPP__

#include "doris_system_info.hpp"
#include "ggdatetime/dtcalendar.hpp"
#include <fstream>
#include <vector>

namespace ids {

/// @class TimeReferenceStation
/// This class holds time reference stations (beacons) recorded in a DORIS
/// RINEX files. These stations are marked with the 'TIME REF STATION' in the
/// RINEX header.
/// @note The m_station_code variable, aka the internal number used in data
///       records, must correspond to a station in the 'STATION REFERENCE'
///       field.
struct TimeReferenceStation {
  /// Internal number used in data records
  char m_station_code[3];
  /// Bias of the time beacon reference VS TAI reference time, unit 1
  /// microsecond
  double m_bias,
      /// Time beacon reference shift unit 1e-14 second/second
      m_shift;
}; // TimeReferenceStation

/// @class RinexDataRecordHeader
/// This class holds fields of a data record header line as in DORIS RINEX
/// files.
/// @see RINEX DORIS 3.0 (Issue 1.7)
struct RinexDataRecordHeader {
  /// reference data record epoch
  dso::datetime<dso::nanoseconds> m_epoch;
  /// Receiver clock offset in seconds (optional)
  double m_clock_offset{RECEIVER_CLOCK_OFFSET_MISSING};
  ///< Number of stations observed in current epoch
  int_fast16_t m_num_stations;
  /// Epoch flag
  int_fast8_t m_flag,
      ///< Receiver clock offset flag, 1 if extrapolated, 0 otherwise
      m_clock_flag;
  /// @brief Apply the recorded (in RINEX) clock offset to time m_epoch
  void apply_clock_offset() noexcept {
    dso::nanoseconds off_nsec {static_cast<dso::nanoseconds::underlying_type>( 
      m_clock_offset * dso::nanoseconds::sec_factor<double>())};
    m_epoch.add_seconds(off_nsec);
  }
}; // RinexDataRecordHeader

struct RinexObservationValue {
  RinexObservationValue(double v, char f1, char f2) noexcept
      : m_value(v), m_flag1(f1), m_flag2(f2){};
  /// The actual value parsed from the corresponding RINEX field
  double m_value;
  char m_flag1, m_flag2;
#ifdef DEBUG
  std::string to_string() const {
    return std::to_string(m_value) + std::string("_f1[") + m_flag1 +
           std::string("]_f2[") + m_flag2 + std::string("]");
  }
#endif
}; // RinexObservationValue

struct BeaconObservations {
  std::vector<RinexObservationValue> m_values;
  char m_beacon_id[3];
  explicit BeaconObservations(int size_hint = 10) noexcept {
    m_values.reserve(size_hint);
  }
#ifdef DEBUG
  std::string to_string() const {
    std::string s(m_beacon_id, 3);
    for (const auto &v : m_values) {
      s += std::string("/") + v.to_string();
    }
    return s;
  }
#endif
}; // BeaconObservations

/// @class DorisObsRinex
/// @brief A class to hold DORIS Observation RINEX files for reading.
/// @see RINEX DORIS 3.0 (Issue 1.7),
///      ftp://ftp.ids-doris.org/pub/ids/data/RINEX_DORIS.pdf
class DorisObsRinex {
public:
  /// Let's not write this more than once.
  typedef std::ifstream::pos_type pos_type;

  /// No header line can have more than 80 chars.
  static constexpr int MAX_HEADER_CHARS{81};

  /// No record line can have more than 3+5*16=83 chars
  static constexpr int MAX_RECORD_CHARS{124};

private:
  /// The name of the file
  std::string m_filename;
  /// The infput (file) stream; open at constructor
  std::ifstream m_stream;
  /// RINEX version
  float m_version;
  /// Satellite name
  char m_satellite_name[60];
  /// COSPAR number
  char m_cospar_number[20];
  /// DORIS chain used (chain1 or chain2), exp. “CHAIN1”
  char m_rec_chain[20];
  /// DORIS instrument type; exp. “DGXX”
  char m_rec_type[20];
  /// The software version used on board DORIS/DIODE, exp. “1.00”
  char m_rec_version[20];
  /// The antenna type is “STAREC”
  char m_antenna_type[20];
  /// The antenna number is “DORIS”
  char m_antenna_number[20];
  /// Position of 2 GHz phase center, in the platform reference frame (Units:
  /// Meters, System: ITRS recommended)
  float m_approx_position[3];
  /// The center of mass of the vehicle (for space borne receivers):
  /// CENTER OF MASS: XYZ, defined at the beginning of the mission.
  float m_center_mass[3];
  /// A vector of ObservationCode contained in the RINEX file
  std::vector<ObservationCode> m_obs_codes;
  /// A vector of scale factors corresponding to m_obs_codes (aka they have
  /// the same size with a one-to-one correspondance
  std::vector<int> m_obs_scale_factors;
  /// Datetime of first observation in RINEX
  dso::datetime<dso::nanoseconds> m_time_of_first_obs;
  /// This date corresponds to the day of the first measurement performed on
  /// the first time reference beacon in the DORIS RINEX product, at
  /// 00h 00mn 00s.
  dso::datetime<dso::nanoseconds> m_time_ref_stat;
  /// Constant shift between the date of the 400MHz phase measurement and the
  /// date of the 2GHz phase measurement in microseconds
  double m_l12_date_offset;
  /// List of stations/beacons recorded in file
  std::vector<BeaconStation> m_stations;
  /// List of time-reference stations in file (also included in m_stations)
  std::vector<TimeReferenceStation> m_ref_stations;
  /// Mark the 'END OF HEADER' field (next line is record line)
  pos_type m_end_of_head;
  /// Record lines for each beacon (in data record blocks)
  int m_lines_per_beacon;

  /// @brief Depending on the number of observables, compute the number of
  /// lines needed to hold a full data record. Each data line can hold up to 5
  /// observable values.
  int lines_per_beacon() const noexcept {
    int obs = m_obs_codes.size();
    return 1 + (!(obs % 5) ? (obs / 5 - 1) : (obs / 5));
  }

  /// @brief read and resolve a RINEX header record.
  int read_header() noexcept;

public:
  /// @brief Constructor from filename
  explicit DorisObsRinex(const char *);

  /// @brief Destructor
  ~DorisObsRinex() noexcept;// = default;

  /// @brief Copy not allowed !
  DorisObsRinex(const DorisObsRinex &) = delete;

  /// @brief Assignment not allowed !
  DorisObsRinex &operator=(const DorisObsRinex &) = delete;

  /// @brief Move Constructor.
  DorisObsRinex(DorisObsRinex &&a) noexcept(
      std::is_nothrow_move_constructible<std::ifstream>::value) = default;

  /// @brief Move assignment operator.
  DorisObsRinex &operator=(DorisObsRinex &&a) noexcept(
      std::is_nothrow_move_assignable<std::ifstream>::value) = default;
  
  /// @brief Given a data record header line, resolve it to a
  ///        RinexDataRecordHeader instance.
  /// @param[in]  line A RINEX data record header line
  /// @param[out] hdr A RinexDataRecordHeader; at output it will hold the info
  ///             resolved from the input line.
  /// @return Anything other than 0 denotes an error; in this case, the hdr
  ///         instance may hold erronuous values and should not be used.
  int resolve_data_epoch(const char *line,
                         RinexDataRecordHeader &hdr) const noexcept;

  /// @brief Read next RINEX data block
  /// @param[in] hdr A RinexDataRecordHeader; the data header record (that
  ///                includes epoch and beacon information) read in the
  ///                start of the data block.
  /// param[out] obsvec A vector of BeaconObservations; for each of the beacons
  ///                recorded in the data block, a new entry is appended in the
  ///                obsvec. For each beacon, the corresponding
  ///                BeaconObservations includes all observation types recorded
  ///                in the RINEX header. For any missing values, the default
  ///                value OBSERVATION_VALUE_MISSING is filled in.
  /// @return Anything other than 0, denotes an error
  int read_data_block(RinexDataRecordHeader &hdr,
                      std::vector<BeaconObservations> &obsvec) noexcept;

  /// @brief Skip next RINEX data block
  /// @param[in] hdr A RinexDataRecordHeader; the data header record (that
  ///                includes epoch and beacon information) read in the
  ///                start of the data block.
  void skip_data_block(const RinexDataRecordHeader &hdr) noexcept;

  int get_doppler_counts() noexcept;

  std::vector<BeaconStation> stations() const noexcept { return m_stations; }

  auto ref_datetime() const noexcept { return m_time_ref_stat;}

  std::ifstream& stream() noexcept { return m_stream; }

  std::vector<ObservationCode> observation_codes() const noexcept { return m_obs_codes;}

#ifdef DEBUG
  void read();
#endif
}; // DorisObsRinex

int fit_relative_frequency_offset(char **rinex_fn, int num_rinex) noexcept;

} // namespace ids

#endif
