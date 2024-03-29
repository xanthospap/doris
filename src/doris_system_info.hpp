#ifndef __IDS_DORIS_SYSTEM_INFO_HPP__
#define __IDS_DORIS_SYSTEM_INFO_HPP__

#include <limits>
#include <stdexcept>
#include <stdint.h>
#include <cmath>
#include "antenna_pcv.hpp"

namespace dso {

/// In DORIS RINEX files, the receiver clock offset may be missing for
/// some/all epochs; this value signifies a missing epoch Receiver clock
/// offset value.
static constexpr double RECEIVER_CLOCK_OFFSET_MISSING =
    std::numeric_limits<double>::min();

/// In DORIS RINEX files, the observation value may be missing for some/all
/// epochs; this value signifies a missing observation value.
static constexpr double OBSERVATION_VALUE_MISSING =
    std::numeric_limits<double>::min();

/// @brief the 2 GHz fundamental DORIS frequency (aka S1)
constexpr double DORIS_FREQ1_MHZ = 2.036250e3; // [MHz]

/// @brief the 400 MHz fundamental DORIS frequency (aka U2)
constexpr double DORIS_FREQ2_MHZ = 401.250e0; // [MHz]

/// @brief the (freq1 / freq2) factor (normally used for iono-free l.
/// combination)
constexpr double GAMMA_FACTOR_SQRT = DORIS_FREQ1_MHZ / DORIS_FREQ2_MHZ;

/// @brief the (freq1 / freq2)^2 factor (normally used for iono-free l.
/// combination)
constexpr double GAMMA_FACTOR = GAMMA_FACTOR_SQRT * GAMMA_FACTOR_SQRT;

/// @brief F0, aka USO frequency in Hz
constexpr double USO_F0 = 5e6;

/// @brief Compute the S1 and U2 (aka 2GHz and 400 MHz) nominal frequencies
///        for a DORIS beacon
/// @param[in] shift_factor The beacon's shift factor (e.g. as extracted from 
///        the 'STATION REFERENCE' field from a DORIS RINEX file)
/// @param[out] s1_freq The S1 (aka 2GHz) nominal frequency [Hz]
/// @param[out] u2_freq The U2 (aka 400MHz) nominal frequency in [Hz]
inline
int beacon_nominal_frequency(int shift_factor, double &s1_freq,
                                       double &u2_freq) noexcept {
  const long two26 = std::pow(2, 26);
  constexpr double fac1 = USO_F0 * 0.75e0;
  const double fac2 =
      (USO_F0 * (87e0 * shift_factor)) / (5e0 * static_cast<double>(two26));
  s1_freq = 543e0 * fac1 + 543e0 * fac2;
  u2_freq = 107e0 * fac1 + 107e0 * fac2;
  return 0;
}

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
/// * 1 to denote the S1 DORIS frequency (on 2 GHz), or
/// * 2 to denote the U2 DORIS frequency (on 400 MHz)
/// @warning m_freq MUST be set to 0, if m_type is one of:
///          * frequency_offset,
///          * ground_pressure,
///          * ground_temperature,
///          * ground_humidity
struct ObservationCode {
  ObservationType m_type;
  int_fast8_t m_freq{0};

  /// @brief Constructor; may throw if the frequency is not valid.
  /// @param[in] type The Observation Type
  /// @param[in] freq The corresponding frequency (if any)
  /// @throw std::runtime_error if the passed-in frequency is not valid.
  ObservationCode(ObservationType type, int_fast8_t freq = 0);

  /// @brief get the ObservationType
  auto type() const noexcept { return m_type; }

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

  /// @brief Format to string (string length=2+null terminating char)
  /// @param[in] buffer A char array of length >= 3 chars
  /// @return a pointer to the (formatted) buffer string
  char* to_str(char *buffer) const noexcept;
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

  /// @return the antenna type (see enum GroundAntennaType)
  GroundAntennaType type() const;

  /// @brief Compute iono-free phase center for this antenna, w.r.t the
  ///        Antenna's reference point. 
  /// @return iono-free phase center w.r.t Antenna's reference point [m] 
  ///        such that:
  ///        Iono-Free PC = Antenna RP + iono_free_phase_center()
  /// @see Lemoine etal, 2016, "Precise orbit determination and station
  ///        position estimation using DORIS RINEX data", section 2.5.7
  Eigen::Matrix<double, 3, 1> iono_free_phase_center() const noexcept {
    // the 2GHz Phase center in [mm]
    const auto r_2Ghz = AntennaOffset<GroundAntennaType::Alcatel, 1>::pco() *
                            (type() == GroundAntennaType::Alcatel) +
                        AntennaOffset<GroundAntennaType::Starec_B, 1>::pco() *
                            (type() == GroundAntennaType::Starec_B ||
                             type() == GroundAntennaType::Starec_C);

    // the 400MHz Phase center in [mm]
    const auto r_400MHz = AntennaOffset<GroundAntennaType::Alcatel, 2>::pco() *
                              (type() == GroundAntennaType::Alcatel) +
                          AntennaOffset<GroundAntennaType::Starec_B, 2>::pco() *
                              (type() == GroundAntennaType::Starec_B ||
                               type() == GroundAntennaType::Starec_C);

    // iono-free phase center in [m] w.r.t Antenna Reference Point
    return (r_2Ghz + (r_2Ghz - r_400MHz) / (GAMMA_FACTOR - 1e0)) * 1e-3;
  }

  /// @brief Concatenate Beacon information to a string
  /// @param[in] buffer A char array of length > 51 chars
  /// @return a pointer to the (formatted) buffer string
  char *to_str(char *buffer) const noexcept;

  /// @brief Set instance's member as resolved from a DORIS RINEX record.
  /// @param[in] line A RINEX 'STATION REFERENCE' record line
  /// @return Anything other than 0 denotes an error.
  /// @see RINEX DORIS 3.0 (Issue 1.7)
  int set_from_rinex_line(const char *line) noexcept;
}; // BeaconStation

} // dso

#endif
