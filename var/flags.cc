#include <cassert>
#include <cstdio>
#include <limits>
#include <cinttypes>
#include <cstdint>
#include <iostream>

enum class Sp3Event : uint_fast16_t {
  bad_abscent_position = 0,
  bad_abscent_clock,
  clock_event,
  clock_prediction,
  maneuver,
  orbit_prediction,
  has_pos_stddev,
  has_clk_stddev,
  bad_abscent_velocity,
  bad_abscent_clock_rate,
  has_vel_stddev,
  has_clk_rate_stdev
}; // Sp3Event

namespace details {
  struct __WrapperFlag__ {
    using uitype = uint_fast16_t;
    uitype bits_{0};
  };//__WrapperFlag__
}// details

details::__WrapperFlag__ operator|(Sp3Event e1, Sp3Event e2) noexcept {
  //std::cout<<"\nSetting bits: "<<static_cast<uint_fast16_t>(e1)<<" and "<<static_cast<uint_fast16_t>(e2);
  details::__WrapperFlag__ wf;
  wf.bits_ |= (1 << static_cast<uint_fast16_t>(e1));
  wf.bits_ |= (1 << static_cast<uint_fast16_t>(e2));
  //std::cout<<"\n\tthat is assign to "<<wf.bits_;
  return wf;
}

details::__WrapperFlag__ operator|(details::__WrapperFlag__ e1, Sp3Event e2) noexcept {
  //std::cout<<"\nSetting bits: "<<static_cast<uint_fast16_t>(e1)<<" and "<<static_cast<uint_fast16_t>(e2);
  details::__WrapperFlag__ wf;
  // wf.bits_ |= (1 << e1.bits_);
  wf.bits_ = e1.bits_;
  wf.bits_ |= (1 << static_cast<uint_fast16_t>(e2));
  //std::cout<<"\n\tthat is assign to "<<wf.bits_;
  return wf;
}

static_assert(std::numeric_limits<uint_fast16_t>::digits >
              static_cast<uint_fast16_t>(Sp3Event::has_clk_rate_stdev));

/// @class Sp3Flag A flag to hold all events recorded in the Sp3 for a field
struct Sp3Flag {
  using uitype = uint_fast16_t;
  /// Initialize unmarked
  uitype bits_{0};
  /// Mark flag with an Sp3Event (aka, set the Sp3Event)
  void set(Sp3Event e) noexcept {
    std::cout<<"\nsetting bit " << static_cast<uitype>(e);
    bits_ |= (1 << static_cast<uitype>(e));
  }
  void set(details::__WrapperFlag__ w) noexcept {
    bits_ = w.bits_;
  }
  /// Un-Mark flag with an Sp3Event (aka, unset the Sp3Event)
  void clear(Sp3Event e) noexcept {
    bits_ &= (~(1 << static_cast<uitype>(e)));
  }
  /// Clear all Sp3Event's and reset flag to empty/clean
  void reset() noexcept { bits_ = 0; }
  /// Trigger, aka check if an Sp3Event is set
  bool is_set(Sp3Event e) const noexcept {
    return ((bits_ >> static_cast<uitype>(e)) & 1);
  }
  /// Check if Sp3Flag is clean (no Sp3Event is set)
  bool is_clean() const noexcept { return !bits_; }
}; // Sp3Flag

int main() {
  Sp3Flag flag;

  assert(flag.is_clean());
  std::cout<<"\nFlag is now: " << flag.bits_;
  flag.set(Sp3Event::bad_abscent_position);
  std::cout<<"\nFlag is now: " << flag.bits_;
  flag.set(Sp3Event::clock_event);
  std::cout<<"\nFlag is now: " << flag.bits_;
  
  flag.reset();
  std::cout<<"\nFlag is now: " << flag.bits_;
  flag.set(Sp3Event::bad_abscent_position|Sp3Event::bad_abscent_clock|Sp3Event::clock_event|Sp3Event::has_clk_rate_stdev|Sp3Event::bad_abscent_clock);
  std::cout<<"\nFlag is now: " << flag.bits_;
  assert(flag.is_set(Sp3Event::bad_abscent_position));
  assert(flag.is_set(Sp3Event::bad_abscent_clock));
  assert(flag.is_set(Sp3Event::clock_event));
  assert(!flag.is_set(Sp3Event::clock_prediction));
  assert(!flag.is_set(Sp3Event::maneuver));
  assert(!flag.is_set(Sp3Event::orbit_prediction));
  assert(!flag.is_set(Sp3Event::has_pos_stddev));
  assert(!flag.is_set(Sp3Event::has_clk_stddev));
  assert(!flag.is_set(Sp3Event::bad_abscent_velocity));
  assert(!flag.is_set(Sp3Event::bad_abscent_clock_rate));
  assert(!flag.is_set(Sp3Event::has_vel_stddev));
  assert(flag.is_set(Sp3Event::has_clk_rate_stdev));

  /*
  flag.set(Sp3Event::bad_abscent_position);
  assert(flag.is_set(Sp3Event::bad_abscent_position));
  assert(!flag.is_set(Sp3Event::bad_abscent_clock));
  assert(!flag.is_set(Sp3Event::clock_event));
  assert(!flag.is_set(Sp3Event::clock_prediction));
  assert(!flag.is_set(Sp3Event::maneuver));
  assert(!flag.is_set(Sp3Event::orbit_prediction));

  flag.reset();
  flag.set(Sp3Event::maneuver | Sp3Event::orbit_prediction);
  assert(!flag.is_set(Sp3Event::bad_abscent_position));
  assert(!flag.is_set(Sp3Event::bad_abscent_clock));
  assert(!flag.is_set(Sp3Event::clock_event));
  assert(!flag.is_set(Sp3Event::clock_prediction));
  assert(flag.is_set(Sp3Event::maneuver));
  assert(flag.is_set(Sp3Event::orbit_prediction));
  assert(!flag.is_set(Sp3Event::has_pos_stddev));
  assert(!flag.is_set(Sp3Event::has_clk_stddev));
  assert(!flag.is_set(Sp3Event::bad_abscent_velocity));
  assert(!flag.is_set(Sp3Event::bad_abscent_clock_rate));
  assert(!flag.is_set(Sp3Event::has_vel_stddev));
  assert(!flag.is_set(Sp3Event::has_clk_rate_stdev));
  */

  std::cout<<"\n";
  return 0;
}
