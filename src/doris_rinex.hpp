#ifndef __IDS_DORIS_RINEX_HPP__
#define __IDS_DORIS_RINEX_HPP__

#include "doris_system_info.hpp"
#include "datetime/dtcalendar.hpp"
#include "filters/models.hpp"
#include <datetime/dtfund.hpp>
#include <fstream>
#include <vector>
#include <algorithm>

namespace dso {

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
  ///< Reference date record epoch (note that this time tag refers to the L1
  /// sampling, for L2 you have to apply the 'L2 / L1 DATE OFFSET')
  dso::datetime<dso::nanoseconds> m_epoch;
  ///< Receiver clock offset in seconds (optional)
  double m_clock_offset{RECEIVER_CLOCK_OFFSET_MISSING};
  ///< Number of stations observed in current epoch
  int_fast16_t m_num_stations;
  ///< Epoch flag
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
  /// @brief Get the beacon (internal) 3-character ID as a **non-null** 
  ///        terminating character array
  /// @warning Cannot stress enough, that is **not a null terminating string**
  const char *id() const noexcept {return m_beacon_id;}
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
/// date of the 2GHz phase measurement in microseconds. Positive if the 
/// measurement of phase 400 MHz is performed after the measurement of phase 
/// 2 GHz
double m_l12_date_offset;
/// Epoch, code, and phase are corrected by applying the realtime-derived 
/// receiver clock offset: 1=yes, 0=no; default: 0=no
bool rcv_clock_offs_appl{false};
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

  const std::vector<BeaconStation> &stations() const noexcept {
    return m_stations;
  }

  auto ref_datetime() const noexcept { return m_time_ref_stat; }

  std::ifstream &stream() noexcept { return m_stream; }

  const std::vector<ObservationCode> &observation_codes() const noexcept {
    return m_obs_codes;
  }

  int get_observation_code_index(ObservationCode t) const noexcept;

  /// @brief Given a beacon 3-char identifier (internal to this RINEX), return 
  ///        the beacon's 4-character station code
  /// @param[in] inid The beacon 3-char identifier (internal to this RINEX). 
  ///        Only the first three chars will be considered, no null-terminated
  ///        string required (just an array of three chars)
  /// @return A pointer to the instance's m_stations[].m_station_id, which is
  ///        a non-null terminating string of size 4.
  ///        If no station is matched, the retuned pointer point to NULL
  /// @warning The character array returned here is a **non-null terminated**
  ///        string. Beware  how you use it. See dso::BeaconStation and its
  ///        member variable m_station_id
  const char *beacon_internal_id2id(const char *_3charid) const noexcept;

  const char *beacon_id2internal_id(const char *_4charid) const noexcept;
  
  /// @brief Given a beacon 3-char identifier (internal to this RINEX), return 
  ///        the corresponding BeaconStation instance (stored in the instance)
  /// @param[in] inid The beacon 3-char identifier (internal to this RINEX). 
  ///        Only the first three chars will be considered, no null-terminated
  ///        string required (just an array of three chars)
  /// @return An iterator to the instance's m_stations, set to the 
  ///        BeaconStation instance with internal (3-char)id equal to the one 
  ///        given. If no station is matched, the retuned iterator is set to
  ///        m_stations.cend().
  std::vector<BeaconStation>::const_iterator
  beacon_internal_id2BeaconStation(const char *inid) const noexcept {
    return std::find_if(m_stations.begin(), m_stations.end(),
                        [&](const BeaconStation &bcn) {
                          return (inid[0] == bcn.m_internal_code[0] &&
                                  inid[1] == bcn.m_internal_code[1] &&
                                  inid[2] == bcn.m_internal_code[2]);
                        });
  }

  auto time_of_first_obs() const noexcept {return m_time_of_first_obs;}

  int print_metadata() const noexcept;

  /// @brief Given a beacon (interbal) 3-char identifier, query the RINEX's
  ///        station list and return the frequency shift factor
  ///        For more information, see DORIS RINEX 3 (1.7), sec. 6.16
  /// @param[in] beaconid A 3-char id of the beacon (as identified within this
  ///        RINEX file, e.g. 'D01')
  /// @param[out] k The shift factor for the queried beacon. To actully
  ///        compute the beacon's nominal frequencies, see 
  ///        ids::beacon_nominal_frequency
  /// @return If a values other than 0 is returned, the beacon was not found
  ///        (in the RINEX) and k should not be used.
  int beacon_shift_factor(const char *beaconid, int &k) const noexcept {
    auto it = std::find_if(stations().cbegin(), stations().cend(),
                           [&](const BeaconStation &b) {
                             return (b.m_internal_code[0] == beaconid[0] &&
                                     b.m_internal_code[1] == beaconid[1] &&
                                     b.m_internal_code[2] == beaconid[2]);
                           });
    if (it != stations().cend()) {
      k = it->m_shift_factor;
      return 0;
    }
    return 1;
  }

  /// @brief Check if receiver clock offsets are applied
  bool receiver_clock_offsets_applied() const noexcept {
    return rcv_clock_offs_appl;
  }

  /// @brief Get the 'L2 / L1 DATE OFFSET' in nanoseconds
  dso::nanoseconds l21_date_offset() const noexcept {
    // note that the value extracted from the RINEX file is in microsecond
    // transform it to nanoseconds (integer)
    const double offset_nanoseconds = m_l12_date_offset * 1e3;
    return dso::nanoseconds(
        static_cast<dso::nanoseconds::underlying_type>(offset_nanoseconds));
  }

#ifdef DEBUG
  void read();
#endif
}; // DorisObsRinex

int fit_relative_frequency_offset(char **rinex_fns, int num_rinex,
                                       double sigma_x = 1e-1,
                                       double sigma_vx = 1e-3,
                                       double sigma_z = 1e1) noexcept;
int fit_relative_frequency_offset(
    const std::vector<const char *> &fns,
    dso::PolynomialModel<dso::datetime<dso::nanoseconds>> &fit,
    bool use_tai = false, int every = 4,
    const dso::datetime<dso::nanoseconds> &start =
        dso::datetime<dso::nanoseconds>::min(),
    const dso::datetime<dso::nanoseconds> &end =
        dso::datetime<dso::nanoseconds>::max()) noexcept;
inline int fit_relative_frequency_offset(
    const char * fns,
    dso::PolynomialModel<dso::datetime<dso::nanoseconds>> &fit,
    bool use_tai = false, int every = 4,
    const dso::datetime<dso::nanoseconds> &start =
        dso::datetime<dso::nanoseconds>::min(),
    const dso::datetime<dso::nanoseconds> &end =
        dso::datetime<dso::nanoseconds>::max()) noexcept
{
  std::vector<const char *> v;
  v.push_back(fns);
  return fit_relative_frequency_offset(v, fit, use_tai, every, start, end);
}

/// @brief An iterator to a DorisObsRinex data.
/// Allows easily iterating of the data block in the RINEX file.
struct RinexDataBlockIterator {
  using value_type = RinexDataRecordHeader;
  using pointer = value_type *;
  using reference = value_type &;

  ///< Current header
  RinexDataRecordHeader cheader;
  ///< Obsrvations in block
  std::vector<BeaconObservations> cblock;
  ///< pointer to the RINEX file
  DorisObsRinex *rnx;

  /// @brief Check if current block contains observations from a given beacon,
  ///        given its internal, 3-char id (RINEX-specific).
  /// @return An iterator to the BeaconObservations instance (in the instances
  ///        clock vector) holding measuremets for the given beacon. If no
  ///        measuremets for the beacon are available (in this block), the 
  ///        iterator in invalid, aka cblock.end()
  std::vector<BeaconObservations>::iterator
  contains_beacon(const char *_3char_id) noexcept {
    return std::find_if(cblock.begin(), cblock.end(),
                        [&](const BeaconObservations &obs) {
                          return obs.m_beacon_id[0] == _3char_id[0] &&
                                 obs.m_beacon_id[1] == _3char_id[1] &&
                                 obs.m_beacon_id[2] == _3char_id[2];
                        });
  }

  /// @brief Constructor
  RinexDataBlockIterator(DorisObsRinex *drnx) noexcept : rnx(drnx) {
    cblock.reserve(5);
  };

  /// @brief Get next data block
  /// Will advance cheader to the next header/epoch and fectch the new
  /// observation set for this epoch. Aka, updates cheader and cblock.
  int next() noexcept;

  dso::datetime<dso::nanoseconds> proper_time() const noexcept {
    return cheader.m_epoch;
  }

  /// @brief Get the L1-reference epoch for the observations, corrected for 
  /// receiver clock offset (if not applied)
  dso::datetime<dso::nanoseconds> corrected_l1_epoch() const noexcept {
    dso::datetime<dso::nanoseconds> t = cheader.m_epoch;
    if (!rnx->receiver_clock_offsets_applied()) [[likely]] {
      // clock offset in nanoseconds (from seconds)
      const double clock_offset =
          static_cast<double>(cheader.m_clock_offset !=
                              RECEIVER_CLOCK_OFFSET_MISSING) *
          cheader.m_clock_offset;
      dso::nanoseconds noff(
          static_cast<dso::nanoseconds::underlying_type>(clock_offset * 1e9));
      t.add_seconds(noff);
    }
    return t;
  }

  /// @brief Get the L1-reference epoch for the observations, corrected for 
  /// receiver clock offset
  /// @param[out] tl2 Reference epoch for the L2 observation, corrected for:
  ///                 1. receiver clock offset, and
  ///                 2. L2/L1 date offset
  dso::datetime<dso::nanoseconds>
  corrected_l1_epoch(dso::datetime<dso::nanoseconds> &tl2) const noexcept;
}; // BlockIterator

} // dso

#endif
