#ifndef __DSO_JASON3_QUATERNION_ROT_HPP__
#define __DSO_JASON3_QUATERNION_ROT_HPP__

#include "eigen3/Eigen/Geometry"
#include "jason3.hpp"
#include <datetime/dtcalendar.hpp>
#include <fstream>
#include <stdexcept>

namespace dso {

struct JasonBodyQuaternion {
  ///< TAI MJD
  dso::TwoPartDate tai_mjd;
  ///< Quaternion
  Eigen::Quaternion<double> quaternion;
}; // JasonBodyQuaternion

/// @ref
/// https://ids-doris.org/documents/BC/ancillary/quaternions/jason1_2_3_quaternion_solar_panel.pdf
struct JasonBodyQuaternionFile {
  std::ifstream fin; ///< the input file stream

  /// @brief JasonBodyQuaternionFile constructor from filename
  /// @param[in] fn filename
  JasonBodyQuaternionFile(const char *fn) noexcept : fin(fn) {}

  /// @brief Get the next record off from the input file
  /// @param record On success, the next time/quaternion from the file
  /// @return Anything other than 0, denotes an error
  int get_next(JasonBodyQuaternion &record) noexcept;

  /// @brief Get the next num_records records off from the input file
  /// @param[in] record On success, the next num_records time/quaternion from
  ///        the file. Should have size at least num_records
  /// @param[in] num_records Number of records (lines) requested
  /// @return Anything other than 0, denotes an error
  int get_next(JasonBodyQuaternion *records, int num_records) noexcept;
};

/// @ref
/// https://ids-doris.org/documents/BC/ancillary/quaternions/jason1_2_3_quaternion_solar_panel.pdf
struct JasonPanelQuaternionFile {
  std::ifstream fin;
  JasonPanelQuaternionFile(const char *fn) noexcept : fin(fn) {}
};

/// @brief Number of JasonBodyQuaternion instances buffered in a
///        JasonQuaternionHunter instance
constexpr const int NumQuaternionsInBuffer = 5;

struct JasonQuaternionHunter {
  JasonBodyQuaternionFile bodyin;
  JasonBodyQuaternion bodyq[NumQuaternionsInBuffer];

  // only for debuging
  void dump_buffered_quaternions() const noexcept {
    for (int i = 0; i < NumQuaternionsInBuffer; i++) {
      printf("\tBuffered Quaternion at: %.2f + %.15f\n", bodyq[i].tai_mjd._big,
             bodyq[i].tai_mjd._small);
    }
  }

  /// @brief  Constructor
  /// @param body_fn Name of the quaternion file, parsed into the instance's
  ///        JasonBodyQuaternionFile member variable
  JasonQuaternionHunter(const char *body_fn);

  /// @brief Move buffered quaternions to the left, aka left shift. This means
  /// that the first quaternion will be lost (replaced by the second) and the
  /// last quaternion will be empty.
  /// After the call, we will have:
  /// bodyq[0] <- bodyq[1]
  /// bodyq[1] <- bodyq[2]
  /// ...
  /// bodyq[NumQuaternionsInBuffer-2] <- bodyq[NumQuaternionsInBuffer-1]
  /// bodyq[NumQuaternionsInBuffer-1]
  void left_shift() noexcept {
    // fuck it, consume memory, no loops
    JasonBodyQuaternion newq[NumQuaternionsInBuffer - 1];
    std::memcpy(newq, bodyq + 1,
                sizeof(JasonBodyQuaternion) * (NumQuaternionsInBuffer - 1));
    std::memcpy(bodyq, newq,
                sizeof(JasonBodyQuaternion) * (NumQuaternionsInBuffer - 1));
    return;
  }

  int find_interval(const dso::TwoPartDate &tai_mjd) const noexcept {
    // printf("Requasting for quaternion at : %.2f = %.15f\n", tai_mjd._big,
    // tai_mjd._small); dump_buffered_quaternions(); start searching from the
    // top, aka from last element
    int qindex = NumQuaternionsInBuffer - 2;
    for (int i = qindex; i >= 0; --i) {
      if ((tai_mjd >= bodyq[i].tai_mjd) && (tai_mjd < bodyq[i + 1].tai_mjd)) {
        return i;
      }
    }
    // tai_mjd is out of bounds, prior to first record in buffer
    if (tai_mjd < bodyq[0].tai_mjd)
      return -1;
    // tai_mjd is out of bounds, after the last record in buffer
    if (tai_mjd >= bodyq[NumQuaternionsInBuffer - 1].tai_mjd)
      return NumQuaternionsInBuffer + 1;
    // we should never reach this point
    return -100;
  }

  /// @brief Go through the input file (if needed), to find a suitable,
  ///        consecutive pair of records such that:
  ///        t>=bodyq[0].t and t<bodyq[1].t
  /// The function operates in a forward manner only (aka goes forward in the
  /// file) and has an effect on where the file streams are placed.
  /// @param[in] t Requested datetime as fractional MJD (TAI)
  /// @return Anything other than 0 signals an error
  int set_at(const dso::TwoPartDate &mjd_tai, int &index) noexcept;

  /// @brief Get the quaternion for a given datetime instance, using the
  ///        SLERP interpolation method.
  /// @param[in]  t Requested datetime as fractional MJD (TAI)
  /// @param[out] q The quaternion at time t
  /// @return Always zero
  int get_at(const dso::TwoPartDate &mjd_tai,
             Eigen::Quaternion<double> &q) noexcept;
}; // JasonQuaternionHunter

} // namespace dso

#endif
