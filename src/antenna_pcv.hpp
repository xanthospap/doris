#ifndef __DORIS_ANTENNA_PCV_HPP__
#define __DORIS_ANTENNA_PCV_HPP__

#include "doris_system_info.hpp"
#include "geodesy/units.hpp"

namespace ids {

/// @class AntennaOffsetTraits A kinda traits class holding Antenna specific
///        information. Needs to be specialized based on:
///        * GroundAntennaType (e.g. Alcatel, Starec, etc)
///        * Frequency (1 or 2)
template <GroundAntennaType T, int Freq> struct AntennaOffsetTraits {};

/// @class Specialization of AntennaOffsetTraits for GroundAntennaType=Alcatel,
///        and Frequency=1
template <> struct AntennaOffsetTraits<GroundAntennaType::Alcatel, 1> {
  /// Eccentricities of the mean antenna phase center relative to the antenna
  /// reference point (ARP). North, east and up component (in millimeters).
  static constexpr double offset[3] = {0e0, 0e0, 510e0};
  /// Phase pattern values in millimeters from 0 to 90 deg, with an increment
  /// of 5 deg.
  static constexpr double plarray[] = {
      0.00e0,  2.05e0,  7.24e0,  9.21e0,  6.71e0,  8.14e0,  11.87e0,
      12.48e0, 12.28e0, 13.67e0, 13.91e0, 13.01e0, 13.01e0, 11.87e0,
      9.70e0,  7.94e0,  4.99e0,  0.41e0,  -3.93}; /* 0-90 step=5*/
  /// @brief Get PCV value in mm, geven the zenith angle using linear
  ///        interpolation
  /// @param[in] zenith Zenith angle in radians
  /// @param[out] out_of_bounds If not 0 (at output), the given zenith angle
  ///            is out of bounds (aka <0 or >90)
  /// @return PCV value for this Antenna/Frequency in mm. If out_of_bounds set,
  ///         the functionwill return 0e0.
  static constexpr double pcv(double zenith, int &out_of_bounds) noexcept {
    const double angle = dso::rad2deg<double>(zenith);
    out_of_bounds = 0;
    if (zenith < 0e0 || zenith > dso::DPI / 2e0) {
      return 0e0;
      out_of_bounds = 1;
    }
    int index = std::floor(angle / 5e0);
#ifdef DEBUG
    assert(index >= 0 && index < static_cast<int>(sizeof(plarray) - 1));
#endif
    return plarray[index] + (angle - 5e0 * index) *
                                ((plarray[index + 1] - plarray[index]) / 5e0);
  }
};

/// @class Specialization of AntennaOffsetTraits for GroundAntennaType=Alcatel,
///        and Frequency=2
template <> struct AntennaOffsetTraits<GroundAntennaType::Alcatel, 2> {
  /// Eccentricities of the mean antenna phase center relative to the antenna
  /// reference point (ARP). North, east and up component (in millimeters).
  static constexpr double offset[3] = {0e0, 0e0, 335e0};
  /// @brief Antenna/Frequency pair has no PCV information; always return 0
  static constexpr double pcv(double zenith, int &out_of_bounds) noexcept {
    out_of_bounds = (zenith < 0e0) + (zenith > dso::DPI / 2e0);
    return 0e0;
  }
};

/// @class Specialization of AntennaOffsetTraits for GroundAntennaType=Starec_B,
///        and Frequency=1
template <> struct AntennaOffsetTraits<GroundAntennaType::Starec_B, 1> {
  /// Eccentricities of the mean antenna phase center relative to the antenna
  /// reference point (ARP). North, east and up component (in millimeters).
  static constexpr double offset[3] = {0e0, 0e0, 487e0};
  /// Phase pattern values in millimeters from 0 to 90 deg, with an increment
  /// of 5 deg.
  static constexpr double plarray[] = {
      0.00e0,  0.06e0,  -0.32e0, -1.12e0, -2.87e0, -4.02e0, -3.44e0,
      -2.15e0, -1.73e0, -1.73e0, -0.08e0, 1.37e0,  2.20e0,  5.37e0,
      7.02e0,  10.70e0, 13.86e0, 17.27e0, 22.37}; /* 0-90 step=5*/
  /// @brief Get PCV value in mm, geven the zenith angle.
  /// @param[in] zenith Zenith angle in radians
  /// @param[out] out_of_bounds If not 0 (at output), the given zenith angle
  ///            is out of bounds (aka <0 or >90)
  /// @return PCV value for this Antenna/Frequency in mm. If out_of_bounds set,
  ///         the functionwill return 0e0.
  static constexpr double pcv(double zenith, int &out_of_bounds) noexcept {
    const double angle(dso::rad2deg<double>(zenith));
    out_of_bounds = 0;
    if (zenith < 0e0 || zenith > dso::DPI / 2e0) {
      out_of_bounds = 1;
      return 0e0;
    }
    int index = std::floor(angle / 5e0);
#ifdef DEBUG
    assert(index >= 0 && index < static_cast<int>(sizeof(plarray) - 1));
#endif
    return plarray[index] + (angle - 5e0 * index) *
                                ((plarray[index + 1] - plarray[index]) / 5e0);
  }
};

/// @class Specialization of AntennaOffsetTraits for GroundAntennaType=Starec_B,
///        and Frequency=2
template <> struct AntennaOffsetTraits<GroundAntennaType::Starec_B, 2> {
  /// Eccentricities of the mean antenna phase center relative to the antenna
  /// reference point (ARP). North, east and up component (in millimeters).
  static double constexpr offset[3] = {0e0, 0e0, 0e0};
  /// @brief Antenna/Frequency pair has no PCV information; always return 0
  static constexpr double pcv(double zenith, int &out_of_bounds) noexcept {
    out_of_bounds = (zenith < 0e0) + (zenith > dso::DPI / 2e0);
    return 0e0;
  }
};

/// @class Specialization of AntennaOffsetTraits for GroundAntennaType=Starec_C,
///        and Frequency=1
template <> struct AntennaOffsetTraits<GroundAntennaType::Starec_C, 1> {
  /// Eccentricities of the mean antenna phase center relative to the antenna
  /// reference point (ARP). North, east and up component (in millimeters).
  static constexpr double offset[3] = {0e0, 0e0, 487e0};
  /// Phase pattern values in millimeters from 0 to 90 deg, with an increment
  /// of 5 deg.
  static constexpr double plarray[] = {
      0.00e0,  0.06e0,  -0.32e0, -1.12e0, -2.87e0, -4.02e0, -3.44e0,
      -2.15e0, -1.73e0, -1.73e0, -0.08e0, 1.37e0,  2.20e0,  5.37e0,
      7.02e0,  10.70e0, 13.86e0, 17.27e0, 22.37}; /* 0-90 step=5*/
  /// @brief Get PCV value in mm, geven the zenith angle.
  /// @param[in] zenith Zenith angle in radians
  /// @param[out] out_of_bounds If not 0 (at output), the given zenith angle
  ///            is out of bounds (aka <0 or >90)
  /// @return PCV value for this Antenna/Frequency in mm. If out_of_bounds set,
  ///         the functionwill return 0e0.
  static constexpr double pcv(double zenith, int &out_of_bounds) noexcept {
    const double angle(dso::rad2deg<double>(zenith));
    out_of_bounds = 0;
    if (zenith < 0e0 || zenith > dso::DPI / 2e0) {
      out_of_bounds = 1;
      return 0e0;
    }
    int index = std::floor(angle / 5e0);
#ifdef DEBUG
    assert(index >= 0 && index < static_cast<int>(sizeof(plarray) - 1));
#endif
    return plarray[index] + (angle - 5e0 * index) *
                                ((plarray[index + 1] - plarray[index]) / 5e0);
  }
};

/// @class Specialization of AntennaOffsetTraits for GroundAntennaType=Starec_C,
///        and Frequency=2
template <> struct AntennaOffsetTraits<GroundAntennaType::Starec_C, 2> {
  /// Eccentricities of the mean antenna phase center relative to the antenna
  /// reference point (ARP). North, east and up component (in millimeters).
  static constexpr double offset[3] = {0e0, 0e0, 0e0};
  /// @brief Antenna/Frequency pair has no PCV information; always return 0
  static constexpr double pcv(double zenith, int &out_of_bounds) noexcept {
    out_of_bounds = (zenith < 0e0) + (zenith > dso::DPI / 2e0);
    return 0e0;
  }
};

/// @class AntennaOffset
/// At construction, use template parameters to use a specific Antenna type
/// and frequency (e.g. AntennaOffset<GroundAntennaType::Alcatel, 1>).
/// Allows interface with the specialized traits structs:
/// template<GroundAntennaType T, int Freq> struct AntennaOffsetTraits
/// You can get offsets (north, east and up) and pcv values (if any).
template <GroundAntennaType T, int Freq>
struct AntennaOffset : AntennaOffsetTraits<T, Freq> {
  /// @brief North eccentricity of the mean antenna phase center in mm
  static constexpr double dnorth() noexcept {
    return AntennaOffsetTraits<T, Freq>::offset[0];
  }
  /// @brief East eccentricity of the mean antenna phase center in mm
  static constexpr double deast() noexcept {
    return AntennaOffsetTraits<T, Freq>::offset[1];
  }
  /// @brief Up eccentricity of the mean antenna phase center in mm
  static constexpr double dup() noexcept {
    return AntennaOffsetTraits<T, Freq>::offset[2];
  }
  /// @brief Get PCV value in mm, geven the zenith angle.
  /// @param[in] zenith Zenith angle in radians
  /// @param[out] out_of_bounds If not 0 (at output), the given zenith angle
  ///            is out of bounds (aka <0 or >90)
  /// @return PCV value for this Antenna/Frequency in mm. If out_of_bounds set,
  ///         the functionwill return 0e0. If the Antenna/Frequency pair has
  ///         no relevant information (aka no PCV values), 0 is returned.
  static constexpr double pcv(double zenithdeg, int &out_of_bounds) noexcept {
    return AntennaOffsetTraits<T, Freq>::pcv(zenithdeg, out_of_bounds);
  }
};

} // namespace ids
#endif
