#include "sp3.hpp"
#include <cassert>
#include <cstdio>

using namespace ids;

int main() {
  Sp3Flag flag;

  assert(flag.is_clean());

  flag.set(Sp3Event::bad_abscent_position);
  assert(flag.is_set(Sp3Event::bad_abscent_position));
  assert(!flag.is_set(Sp3Event::bad_abscent_clock));
  assert(!flag.is_set(Sp3Event::clock_event));
  assert(!flag.is_set(Sp3Event::clock_prediction));
  assert(!flag.is_set(Sp3Event::maneuver));
  assert(!flag.is_set(Sp3Event::orbit_prediction));
  
  flag.set(Sp3Event::bad_abscent_clock);
  assert(flag.is_set(Sp3Event::bad_abscent_position));
  assert(flag.is_set(Sp3Event::bad_abscent_clock));
  assert(!flag.is_set(Sp3Event::clock_event));
  assert(!flag.is_set(Sp3Event::clock_prediction));
  assert(!flag.is_set(Sp3Event::maneuver));
  assert(!flag.is_set(Sp3Event::orbit_prediction));
  
  flag.set(Sp3Event::clock_event);
  assert(flag.is_set(Sp3Event::bad_abscent_position));
  assert(flag.is_set(Sp3Event::bad_abscent_clock));
  assert(flag.is_set(Sp3Event::clock_event));
  assert(!flag.is_set(Sp3Event::clock_prediction));
  assert(!flag.is_set(Sp3Event::maneuver));
  assert(!flag.is_set(Sp3Event::orbit_prediction));
  
  flag.set(Sp3Event::clock_prediction);
  assert(flag.is_set(Sp3Event::bad_abscent_position));
  assert(flag.is_set(Sp3Event::bad_abscent_clock));
  assert(flag.is_set(Sp3Event::clock_event));
  assert(flag.is_set(Sp3Event::clock_prediction));
  assert(!flag.is_set(Sp3Event::maneuver));
  assert(!flag.is_set(Sp3Event::orbit_prediction));
  
  flag.set(Sp3Event::maneuver);
  assert(flag.is_set(Sp3Event::bad_abscent_position));
  assert(flag.is_set(Sp3Event::bad_abscent_clock));
  assert(flag.is_set(Sp3Event::clock_event));
  assert(flag.is_set(Sp3Event::clock_prediction));
  assert(flag.is_set(Sp3Event::maneuver));
  assert(!flag.is_set(Sp3Event::orbit_prediction));
  
  flag.set(Sp3Event::orbit_prediction);
  assert(flag.is_set(Sp3Event::bad_abscent_position));
  assert(flag.is_set(Sp3Event::bad_abscent_clock));
  assert(flag.is_set(Sp3Event::clock_event));
  assert(flag.is_set(Sp3Event::clock_prediction));
  assert(flag.is_set(Sp3Event::maneuver));
  assert(flag.is_set(Sp3Event::orbit_prediction));

  flag.clear(Sp3Event::orbit_prediction);
  assert(flag.is_set(Sp3Event::bad_abscent_position));
  assert(flag.is_set(Sp3Event::bad_abscent_clock));
  assert(flag.is_set(Sp3Event::clock_event));
  assert(flag.is_set(Sp3Event::clock_prediction));
  assert(flag.is_set(Sp3Event::maneuver));
  assert(!flag.is_set(Sp3Event::orbit_prediction));

  flag.clear(Sp3Event::bad_abscent_position);
  assert(!flag.is_set(Sp3Event::bad_abscent_position));
  assert(flag.is_set(Sp3Event::bad_abscent_clock));
  assert(flag.is_set(Sp3Event::clock_event));
  assert(flag.is_set(Sp3Event::clock_prediction));
  assert(flag.is_set(Sp3Event::maneuver));
  assert(!flag.is_set(Sp3Event::orbit_prediction));
  
  flag.clear(Sp3Event::clock_prediction);
  assert(!flag.is_set(Sp3Event::bad_abscent_position));
  assert(flag.is_set(Sp3Event::bad_abscent_clock));
  assert(flag.is_set(Sp3Event::clock_event));
  assert(!flag.is_set(Sp3Event::clock_prediction));
  assert(flag.is_set(Sp3Event::maneuver));
  assert(!flag.is_set(Sp3Event::orbit_prediction));
  
  flag.clear(Sp3Event::bad_abscent_clock);
  assert(!flag.is_set(Sp3Event::bad_abscent_position));
  assert(!flag.is_set(Sp3Event::bad_abscent_clock));
  assert(flag.is_set(Sp3Event::clock_event));
  assert(!flag.is_set(Sp3Event::clock_prediction));
  assert(flag.is_set(Sp3Event::maneuver));
  assert(!flag.is_set(Sp3Event::orbit_prediction));
  
  flag.reset();
  assert(!flag.is_set(Sp3Event::bad_abscent_position));
  assert(!flag.is_set(Sp3Event::bad_abscent_clock));
  assert(!flag.is_set(Sp3Event::clock_event));
  assert(!flag.is_set(Sp3Event::clock_prediction));
  assert(!flag.is_set(Sp3Event::maneuver));
  assert(!flag.is_set(Sp3Event::orbit_prediction));

  assert(flag.is_clean());

  flag.set(Sp3Event::bad_abscent_position |
    Sp3Event::bad_abscent_clock |
    Sp3Event::clock_event |
    Sp3Event::clock_prediction |
    Sp3Event::maneuver |
    Sp3Event::orbit_prediction |
    Sp3Event::bad_abscent_position); // duplicate !
  assert(flag.is_set(Sp3Event::bad_abscent_position));
  assert(flag.is_set(Sp3Event::bad_abscent_clock));
  assert(flag.is_set(Sp3Event::clock_event));
  assert(flag.is_set(Sp3Event::clock_prediction));
  assert(flag.is_set(Sp3Event::maneuver));
  assert(flag.is_set(Sp3Event::orbit_prediction));
  assert(!flag.is_set(Sp3Event::has_pos_stddev));
  assert(!flag.is_set(Sp3Event::has_clk_stddev));
  assert(!flag.is_set(Sp3Event::bad_abscent_velocity));
  assert(!flag.is_set(Sp3Event::bad_abscent_clock_rate));
  assert(!flag.is_set(Sp3Event::has_vel_stddev));
  assert(!flag.is_set(Sp3Event::has_clk_rate_stdev));

  printf("All ok!\n");
  return 0;
}
