#ifndef __IDS_DORIS_SYSTEM_INFO_HPP__
#define __IDS_DORIS_SYSTEM_INFO_HPP__

#include <limits>
#include <stdexcept>
#include <stdint.h>

namespace ids {

/// @enum ObservationType
/// DORIS Observation Types as defined in RINEX DORIS 3.0 (Issue 1.7)
enum class ObservationType : int_fast8_t {
  phase,              ///< L
  pseudorange,        ///< C
  power_level,        ///< W power level received at each frequency, unit dBm
  frequency_offset,   ///< F relative frequency offset of the receiver’s
                      ///< oscillator (f-f0) / f0, unit 10e-11
  ground_pressure,    ///< P ground pressure at the station, unit 100 Pa (mBar)
  ground_temperature, ///< T ground temperature at the station, unit degree
                      ///< Celsius
  ground_humidity,    ///< H ground humidity at the station, unit percent
};                    // ObservationType

/// @brief Translate an ObservationType to a character.
/// @param[in] obs An ObservationType to translate to char
/// @throw std::runtime_error if no corresponding char is found for the given
///        ObservationType.
char ObservationType_to_char(ObservationType obs);

/// @brief Translate a character to an ObservationType.
/// @param[in] c A char to translate to an ObservationType
/// @throw std::runtime_error if no corresponding ObservationType is found to
///        match the given character.
ObservationType char_to_observationType(char c);

/// @brief Check if a ObservationType needs corresponding frequency information
/// @param[in] type An ObservationType
/// @return True if type is any of phase, pseudorange or power_level; false
///         otherwise
bool observationType_has_frequency(ObservationType type) noexcept;

/// @brief Observation Code as defined in RINEX DORIS 3.0 (Issue 1.7)
/// An Observation code is actually a colletion of an Observation Type and
/// (if needed) a frequency. Frequency numbers are only relevant for
/// Observation Types of type ObservationType::phase,
/// ObservationType::pseudorange and ObservationType::power_level. In any other
/// case, m_freq is irrelevant and set to 0.
/// Frequency is defined by an integer, which can be:
/// * 1 to denote the S1 DORIS frequency, or
/// * 2 to denote the U2 DORIS frequency
/// @warning m_freq MUST be set to 0, if m_type is one of:
///          * frequency_offset,
///          * ground_pressure,
///          * ground_temperature,
///          * ground_humidity
class ObservationCode {
private:
  ObservationType m_type;
  int_fast8_t m_freq{0};

public:
  /// @brief Constructor; may throw if the frequency is not valid.
  /// @param[in] type The Observation Type
  /// @param[in] freq The corresponding frequency (if any)
  /// @throw std::runtime_error if the passed-in frequency is not valid.
  ObservationCode(ObservationType type, int_fast8_t freq = 0);

  /// @brief Equality comparisson; checks both ObservationType and frequency.
  bool operator==(const ObservationCode &oc) const noexcept {
    return m_type == oc.m_type && m_freq == oc.m_freq;
  }

  /// @brief InEquality comparisson; checks both ObservationType and frequency.
  bool operator!=(const ObservationCode &oc) const noexcept {
    return !(this->operator==(oc));
  }

  /// @brief Check if the ObservationCode has a corresponding frequency
  bool has_frequency() const noexcept;
}; // ObservationCode

/// @brief A station (aka beacon) as defined in RINEX DORIS 3.0 (Issue 1.7)
struct BeaconStation {
  /// Internal number used in data records
  char m_internal_code[3];
  /// 4-character station code
  char m_station_id[4];
  /// Station name
  char m_station_name[30];
  /// DOMES number
  char m_station_domes[10];
  /// Type, 1 for beacon 1.0, 2 for beacon 2.0 or 3 for beacon 3.0
  int_fast8_t m_type;
  /// Frequency shift factor K (signed)
  int m_shift_factor;

  /// @brief Set instance's member as resolved from a DORIS RINEX record.
  /// @param[in] line A RINEX 'STATION REFERENCE' record line
  /// @return Anything other than 0 denotes an error.
  /// @see RINEX DORIS 3.0 (Issue 1.7)
  int set_from_rinex_line(const char *line) noexcept;
}; // BeaconStation

} // namespace ids

#endif
