#ifndef __DSO_JASON3_QUATERNION_ROT_HPP__
#define __DSO_JASON3_QUATERNION_ROT_HPP__

#include "base_error.hpp"
#include "eigen3/Eigen/Geometry"
#include "jason3.hpp"
#include <fstream>
#include <stdexcept>

namespace dso {

struct JasonBodyQuaternion {
  /* TAI MJD (reference time for quaternion in TAI) */
  dso::TwoPartDate tai_mjd;
  /* the attitude Quaternion */
  Eigen::Quaternion<double> quaternion;
}; /* JasonBodyQuaternion */

struct JasonSolarArrayRotation {
  /* TAI MJD (reference time in TAI) */
  dso::TwoPartDate tai_mjd;
  /* Angular position of the left Solar Array [rad]; Rotation axis = -Y*/
  double left_array;
  /* Angular position of the right Solar Array [rad]; Rotation axis = +Y */
  double right_array;
}; /* JasonSolarArrayRotation */

/* @ref
 * https://ids-doris.org/documents/BC/ancillary/quaternions/jason1_2_3_quaternion_solar_panel.pdf
 */
struct JasonBodyQuaternionFile {
  /* Body quaternion file stream */
  std::ifstream fin;

  /* @brief JasonBodyQuaternionFile constructor from filename
   * @param[in] fn filename
   */
  JasonBodyQuaternionFile(const char *fn) noexcept : fin(fn) {}

  /* @brief Get the next record off from the input file
   * @param record On success, the next time/quaternion from the file
   * @return Anything other than 0, denotes an error
   */
  dso::iStatus get_next(JasonBodyQuaternion &record) noexcept;

  /* @brief Get the next num_records records off from the input file
   * @param[in] record On success, the next num_records time/quaternion from
   *        the file. Should have size at least num_records
   * @param[in] num_records Number of records (lines) requested
   * @return Anything other than 0, denotes an error
   */
  dso::iStatus get_next(JasonBodyQuaternion *records, int num_records) noexcept;
}; /* JasonBodyQuaternionFile */

/* @ref
 * https://ids-doris.org/documents/BC/ancillary/quaternions/jason1_2_3_quaternion_solar_panel.pdf
 */
struct JasonSolarArrayFile {
  /* Solar array file stream */
  std::ifstream fin;

  JasonSolarArrayFile(const char *fn) noexcept : fin(fn) {}

  /* @brief Get the next record off from the input file
   * @param record On success, the next time/quaternion from the file
   * @return Anything other than 0, denotes an error
   */
  dso::iStatus get_next(JasonSolarArrayRotation &record) noexcept;

  /* @brief Get the next num_records records off from the input file
   * @param[in] record On success, the next num_records time/array angles from
   *        the file. Should have size at least num_records
   * @param[in] num_records Number of records (lines) requested
   * @return Anything other than 0, denotes an error
   */
  dso::iStatus get_next(JasonSolarArrayRotation *records,
                        int num_records) noexcept;
}; /* JasonSolarArrayFile */

/// @ref
/// https://ids-doris.org/documents/BC/ancillary/quaternions/jason1_2_3_quaternion_solar_panel.pdf
struct JasonPanelQuaternionFile {
  /* Panel quaternion file stream */
  std::ifstream fin;
  /* @brief JasonPanelQuaternionFile constructor from filename
   * @param[in] fn filename
   */
  JasonPanelQuaternionFile(const char *fn) noexcept : fin(fn) {}
}; /* JasonPanelQuaternionFile */

/* @brief Number of JasonBodyQuaternion instances buffered in a
 *        JasonQuaternionHunter instance
 */
constexpr const int NumQuaternionsInBuffer = 20;

/* @brief Attitude of Jason satellites using CNES-distributed quaternion files
 *        (body frame part)
 */
class JasonQuaternionHunter {
  JasonBodyQuaternionFile bodyin;
  JasonBodyQuaternion bodyq[NumQuaternionsInBuffer];

  /* @brief Keep on reading and buffering quaternions untill we find a
   * matching time interval
   */
  dso::iStatus read_untill_buffered(const dso::TwoPartDate &tai_mjd,
                                    int &index) noexcept;

  /* @brief Move buffered quaternions to the left, aka left shift. This means
   * that the first quaternion will be lost (replaced by the second) and the
   * last quaternion will be empty.
   * After the call, we will have:
   * bodyq[0] <- bodyq[1]
   * bodyq[1] <- bodyq[2]
   * ...
   * bodyq[NumQuaternionsInBuffer-2] <- bodyq[NumQuaternionsInBuffer-1]
   * bodyq[NumQuaternionsInBuffer-1]
   */
  void left_shift(int pos = 1) noexcept {
    /* fuck it, consume memory, no loops */
    JasonBodyQuaternion newq[NumQuaternionsInBuffer /* - pos*/];
    std::memcpy(newq, bodyq + pos,
                sizeof(JasonBodyQuaternion) * (NumQuaternionsInBuffer - pos));
    std::memcpy(bodyq, newq,
                sizeof(JasonBodyQuaternion) * (NumQuaternionsInBuffer - pos));
    return;
  }

  /* @brief Search through the buffered quaternions to find a suitable
   * interval for the given epoch
   */
  int find_interval(const dso::TwoPartDate &tai_mjd) const noexcept;

  /* @brief Go through the input file (if needed), to find a suitable,
   *        consecutive pair of records such that:
   *        t>=bodyq[0].t and t<bodyq[1].t
   * The function operates in a forward manner only (aka goes forward in the
   * file) and has an effect on where the file streams are placed.
   * @param[in] t Requested datetime as fractional MJD (TAI)
   * @return Anything other than 0 signals an error
   */
  dso::iStatus set_at(const dso::TwoPartDate &mjd_tai, int &index) noexcept;

public:
  // only for debuging
  void dump_buffered_quaternions() const noexcept {
    fprintf(stderr, "--->\n");
    for (int i = 0; i < NumQuaternionsInBuffer; i++) {
      fprintf(stderr, "\tBuffered Quaternion at: %.2f + %.15f\n",
              bodyq[i].tai_mjd.big(), bodyq[i].tai_mjd.small());
    }
    fprintf(stderr, "--->\n");
  }

  /* @brief  Constructor
   * @param body_fn Name of the quaternion file, parsed into the instance's
   *        JasonBodyQuaternionFile member variable
   */
  JasonQuaternionHunter(const char *body_fn);

  /* @brief Get the quaternion for a given datetime instance, using the
   *        SLERP interpolation method.
   * @param[in]  t Requested datetime as fractional MJD (TAI)
   * @param[out] q The quaternion at time t
   */
  dso::iStatus get_at(const dso::TwoPartDate &mjd_tai,
                      Eigen::Quaternion<double> &q) noexcept;
}; /* JasonQuaternionHunter */

/* @brief Attitude of Jason satellites using CNES-distributed quaternion files
 *        (solar array part)
 */
class JasonSolarArrayHunter {
  JasonSolarArrayFile bodyin;
  JasonSolarArrayRotation rots[NumQuaternionsInBuffer];

  /* @brief Keep on reading and buffering solar angles untill we find a
   * matching time interval
   */
  dso::iStatus read_untill_buffered(const dso::TwoPartDate &tai_mjd,
                                    int &index) noexcept;

  /* @brief Move buffered angles to the left, aka left shift. This means
   * that the first angles will be lost (replaced by the second) and the
   * last angles will be empty.
   * After the call, we will have:
   * rots[0] <- rots[1]
   * rots[1] <- rots[2]
   * ...
   * rots[NumQuaternionsInBuffer-2] <- rots[NumQuaternionsInBuffer-1]
   * rots[NumQuaternionsInBuffer-1]
   */
  void left_shift(int pos = 1) noexcept {
    /* fuck it, consume memory, no loops */
    JasonSolarArrayRotation newq[NumQuaternionsInBuffer];
    std::memcpy(newq, rots + pos,
                sizeof(JasonSolarArrayRotation) *
                    (NumQuaternionsInBuffer - pos));
    std::memcpy(rots, newq,
                sizeof(JasonSolarArrayRotation) *
                    (NumQuaternionsInBuffer - pos));
    return;
  }

  /* @brief Search through the buffered solar angles to find a suitable
   * interval for the given epoch
   */
  int find_interval(const dso::TwoPartDate &tai_mjd) const noexcept;

  /* @brief Go through the input file (if needed), to find a suitable,
   *        consecutive pair of records such that:
   *        t>=rots[0].t and t<rots[1].t
   * The function operates in a forward manner only (aka goes forward in the
   * file) and has an effect on where the file streams are placed.
   * @param[in] t Requested datetime as fractional MJD (TAI)
   * @return Anything other than 0 signals an error
   */
  dso::iStatus set_at(const dso::TwoPartDate &mjd_tai, int &index) noexcept;

public:
  /* @brief Constructor
   * @param array_fn Name of the solar array file
   */
  JasonSolarArrayHunter(const char *array_fn);

  /* @brief Get the array angles for a given datetime instance, using the
   *        linear interpolation method.
   * @param[in]  t Requested datetime as fractional MJD (TAI)
   */
  dso::iStatus get_at(const dso::TwoPartDate &mjd_tai, double &left_angle,
                      double &right_angle) noexcept;
}; /* JasonSolarArrayHunter */

} /* namespace dso */

#endif
