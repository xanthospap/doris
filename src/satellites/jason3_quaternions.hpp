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

struct JasonSolarArrayAngles {
  ///< TAI MJD
  dso::TwoPartDate tai_mjd;
  ///< Left and Right panel agles
  /* See description in file: 
   * https://ids-doris.org/documents/BC/ancillary/quaternions/jason1_2_3_quaternion_solar_panel.pdf
   * Left Solar Array:
   *  Rotation axis = -Y; Angle=0 when solar array's normal direction is oriented
   * along -X satellite's axis. Angle is counted positive for
   * counter clockwise rotation (right hand rule)
   * Right Solar Array:
   *  Rotation axis = +Y; Angle=0 when solar array's normal direction is oriented
   * along -X satellite's axis. Angle is counted positive for
   * counter clockwise rotation (right hand rule)
   */
  double left,right;
};// JasonSolarArrayAngles

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
struct JasonSolarArrayAnglesFile {
  std::ifstream fin;
  JasonSolarArrayAnglesFile(const char *fn) noexcept : fin(fn) {}
  int get_next(JasonSolarArrayAngles *records, int num_records) noexcept;
  int get_next(JasonSolarArrayAngles &record) noexcept;
};

/// @brief Number of JasonBodyQuaternion instances buffered in a
///        JasonQuaternionHunter instance
constexpr const int NumQuaternionsInBuffer = 20;

struct JasonSolarArrayHunter {
  JasonSolarArrayAnglesFile bodyin;
  JasonSolarArrayAngles angles[NumQuaternionsInBuffer];
  
  /// @brief  Constructor
  /// @param array_fn Name of the solar array file, parsed into the instance's
  JasonSolarArrayHunter(const char *body_fn);
  int find_interval(const dso::TwoPartDate &tai_mjd) const noexcept;
  void left_shift() noexcept {
    // fuck it, consume memory, no loops
    JasonBodyQuaternion newq[NumQuaternionsInBuffer - 1];
    std::memcpy(newq, angles + 1,
                sizeof(JasonSolarArrayAngles) * (NumQuaternionsInBuffer - 1));
    std::memcpy(angles, newq,
                sizeof(JasonSolarArrayAngles) * (NumQuaternionsInBuffer - 1));
    return;
  }
  int set_at(const dso::TwoPartDate &mjd_tai, int &index) noexcept;
  int get_at(const dso::TwoPartDate &mjd_tai,
             double &lefty, double &righty) noexcept;
};// JasonSolarArrayHunter

struct JasonQuaternionHunter {
  JasonBodyQuaternionFile bodyin;
  JasonBodyQuaternion bodyq[NumQuaternionsInBuffer];

  // only for debuging
  void dump_buffered_quaternions() const noexcept {
    fprintf(stderr, "--->\n");
    for (int i = 0; i < NumQuaternionsInBuffer; i++) {
      fprintf(stderr, "\tBuffered Quaternion at: %.2f + %.15f\n",
              bodyq[i].tai_mjd._big, bodyq[i].tai_mjd._small);
    }
    fprintf(stderr, "--->\n");
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

  int find_interval(const dso::TwoPartDate &tai_mjd) const noexcept;

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

struct JasonAttitude {
  JasonQuaternionHunter q;
  JasonSolarArrayHunter s;
  
  JasonAttitude(const char *fnbody, const char *fnarray) noexcept
      : q(fnbody), s(fnarray){};
  
  int get_at(const dso::TwoPartDate &mjd_tai, Eigen::Quaternion<double> &qtrn,
             double &left, double &right) noexcept {
    int error = q.get_at(mjd_tai, qtrn);
    error += s.get_at(mjd_tai, left, right);
    return error;
  }
  int get_at(const dso::TwoPartDate &mjd_tai, Eigen::Quaternion<double> &qtrn) noexcept {
    int error = q.get_at(mjd_tai, qtrn);
    return error;
  }
};// JasonAttitude

} // namespace dso

#endif
