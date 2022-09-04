#ifndef __DSO_JASON3_QUATERNION_ROT_HPP__
#define __DSO_JASON3_QUATERNION_ROT_HPP__

#include "jason3.hpp"
#include <fstream>
#include <stdexcept>
#include "eigen3/Eigen/Geometry"

namespace dso {

struct JasonBodyQuaternion {
    ///< datetime, UTC
    dso::datetime<dso::nanoseconds> t{dso::datetime<dso::nanoseconds>::min()};
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
    int get_next(dso::JasonBodyQuaternion *records, int num_records) noexcept;
};

/// @ref
/// https://ids-doris.org/documents/BC/ancillary/quaternions/jason1_2_3_quaternion_solar_panel.pdf
struct JasonPanelQuaternionFile {
    std::ifstream fin;
    JasonPanelQuaternionFile(const char *fn) noexcept : fin(fn) {}
};

struct JasonQuaternionHunter {
  JasonBodyQuaternionFile bodyin;
  JasonBodyQuaternion bodyq[2];

  JasonQuaternionHunter(const char *body_fn);

  /// @brief Go through the input file (if needed), to find a suitable, 
  ///        consecutive pair of records such that:
  ///        t>=bodyq[0].t and t<bodyq[1].t
  /// The function operates in a forward manner only (aka goes forward in the
  /// file) and has an effect on where the file streams are placed.
  /// @param[in] t Requested datetime
  /// @return Anythiong other than 0 signals an error
  int set_at(const dso::datetime<dso::nanoseconds> &t) noexcept;

  /// @brief Get the quaternion for a given datetime instance, using the
  ///        SLERP interpolation method.
  /// @param[in]  t Requested datetime
  /// @param[out] q The quaternion at time t
  /// @return Always zero
  int get_at(const dso::datetime<dso::nanoseconds> &t,
             Eigen::Quaternion<double> &q) noexcept;
}; // JasonQuaternionHunter

} // dso

#endif
