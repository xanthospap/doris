#ifndef __SP3C_IGS_FILE__
#define __SP3C_IGS_FILE__

#include "ggdatetime/dtcalendar.hpp"
#include <fstream>
#ifdef DEBUG
#include "ggdatetime/datetime_write.hpp"
#endif

namespace ids {

enum class Sp3Event : unsigned int {
  bad_abscent_position = 0,
  bad_abscent_clock,
  clock_event,
  clock_prediction,
  maneuver,
  orbit_prediction
};// Sp3Event

static_assert(std::numeric_limits<unsigned char>::digits >
              static_cast<unsigned int>(Sp3Event::orbit_prediction));

struct Sp3Flag {
  unsigned char bits_{0};
  void set(Sp3Event e) noexcept {
    bits_ |= (1 << static_cast<unsigned char>(e));
  }
  void clear(Sp3Event e) noexcept {
    bits_ &= (~(1 << static_cast<unsigned char>(e)));
  }
  void reset() noexcept { bits_ = 0; }
  bool is_set(Sp3Event e) const noexcept {
    return ((bits_ >> static_cast<unsigned char>(e)) & 1);
  }
  bool is_clean() const noexcept {
    return !bits_;
  }
};// Sp3Flag

class Sp3c {
public:
  /// Let's not write this more than once.
  typedef std::ifstream::pos_type pos_type;

  /// @brief Constructor from filename
  explicit Sp3c(const char *);

  /// @brief Destructor (closing the file is not mandatory, but nevertheless)
  ~Sp3c() noexcept {
    if (__istream.is_open())
      __istream.close();
  }

  /// @brief Copy not allowed !
  Sp3c(const Sp3c &) = delete;

  /// @brief Assignment not allowed !
  Sp3c &operator=(const Sp3c &) = delete;

  /// @brief Move Constructor.
  Sp3c(Sp3c &&a) noexcept(
      std::is_nothrow_move_constructible<std::ifstream>::value) = default;

  /// @brief Move assignment operator.
  Sp3c &operator=(Sp3c &&a) noexcept(
      std::is_nothrow_move_assignable<std::ifstream>::value) = default;

  auto interval() const noexcept { return interval__; }

  auto num_sats() const noexcept { return num_sats__; }

private:
  /// @brief Read sp3c header; assign info
  int read_header() noexcept;

  std::string __filename;  ///< The name of the file
  std::ifstream __istream; ///< The infput (file) stream
  char version__;          ///< the version 'c' or 'd'
  ngpt::datetime<ngpt::microseconds> start_epoch__; ///< Start epoch
  int num_epochs__,              ///< Number of epochs in file
      num_sats__;                ///< Number od SVs in file
  char crd_sys__[6],///< Coordinate system (last char always '\0')
    orb_type__[4], ///< Orbit type (last char always '\0')
    agency__[5],  ///< Agency (last char always '\0')
    time_sys__[4]; ///< Time system (last char always '\0')
  ngpt::microseconds interval__; ///< Epoch interval
  //SATELLITE_SYSTEM __satsys;     ///< satellite system
  pos_type __end_of_head;        ///< Mark the 'END OF HEADER' field
};// Sp3c

}// ids

#endif
