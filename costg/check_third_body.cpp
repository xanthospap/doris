#include "astrodynamics.hpp"
#include "costg_utils.hpp"
#include "egravity.hpp"
#include "orbit_integration.hpp"
#include "planetpos.hpp"
#include <cassert>

constexpr const double GM_Sun = 1.32712442076e20;
constexpr const double GM_Moon = 0.49028010560e13;
[[maybe_unused]] constexpr const double GM_MERCURY = GM_Sun / 6023600.0;
[[maybe_unused]] constexpr const double GM_VENUS = GM_Sun / 408523.71;
[[maybe_unused]] constexpr const double GM_MARS = GM_Sun / 3098708.0;
[[maybe_unused]] constexpr const double GM_JUPITER = GM_Sun / 1047.3486;
[[maybe_unused]] constexpr const double GM_SATURN = GM_Sun / 3497.898;

int main(int argc, char *argv[]) {
  if (argc != 6) {
    fprintf(stderr,
            "Usage: %s [03directTideMoon_icrf.txt] [03directTideSun_icrf.txt]"
            "[00orbit_icrf.txt] [de421.bsp] [naif0012.tls]\n",
            argv[0]);
    return 1;
  }

  /* Load CSPICE/NAIF Kernels */
  dso::cspice::load_if_unloaded_spk(argv[4]);
  dso::cspice::load_if_unloaded_spk(argv[5]);

  /* parse COST-G (acceleration) results, for Moon ... */
  std::vector<costg::CostgAcc> avmoon;
  if (costg::parse_gravity_field(argv[1], avmoon)) {
    fprintf(stderr, "ERROR Failed to read input acceleration file %s\n",
            argv[1]);
    return 1;
  }
  /* ... for Sun */
  std::vector<costg::CostgAcc> avsun;
  if (costg::parse_gravity_field(argv[2], avsun)) {
    fprintf(stderr, "ERROR Failed to read input acceleration file %s\n",
            argv[2]);
    return 1;
  }

  /* parse COST-G state, ICRF */
  std::vector<costg::CostgExtState> vstate;
  if (costg::parse_satellite_state(argv[3], vstate)) {
    fprintf(stderr, "ERROR Failed to read input state file %s\n", argv[3]);
    return 1;
  }

  /* holds position of Sun and Moon (unused) */
  Eigen::Matrix<double, 3, 1> rsun, rmon;

  int error = 0;

  /* loop through reference results ... */
  auto state = vstate.begin();
  auto amoon = avmoon.begin();
  auto asun = avsun.begin();
  while (state != vstate.end()) {
    amoon = std::find_if(amoon, avmoon.end(), [&](const costg::CostgAcc &g) {
      return g.gpst == state->gpst;
    });
    asun = std::find_if(asun, avsun.end(), [&](const costg::CostgAcc &g) {
      return g.gpst == state->gpst;
    });
    if (amoon == avmoon.end() || asun == avsun.end()) {
      fprintf(stderr, "ERROR Failed to find matching gravity result!\n");
      return 1;
    }

    /* gps time to TAI */
    const auto tai = costg::gps2tai(state->gpst);


    /* position of Sun and Moon at current epoch, ICRF */
    if (dso::planet_pos(dso::Planet::SUN, tai.tai2tt(), rsun)) ++error;
    if (dso::planet_pos(dso::Planet::MOON, tai.tai2tt(), rmon)) ++error;

    if (error) return 1;

    /* point-mass acceleration, Moon */
    auto moon_acc = dso::point_mass_accel(GM_Moon, state->r(), rmon);
    
    /* point-mass acceleration, Sun */
    auto sun_acc = dso::point_mass_accel(GM_Sun, state->r(), rsun);

    /* print results */
    printf("[MOON] %.12e %.20e %.20e %.20e %.20e %.20e %.20e\n",
           amoon->gpst.mjd(), amoon->ax, amoon->ay, amoon->az, moon_acc(0),
           moon_acc(1), moon_acc(2));

    printf("[SUN] %.12e %.20e %.20e %.20e %.20e %.20e %.20e\n",
           asun->gpst.mjd(), asun->ax, asun->ay, asun->az, sun_acc(0),
           sun_acc(1), sun_acc(2));

    /* augment state to next */
    ++state;
  }

  return 0;
}
