#include "costg_utils.hpp"
#include "iers2010/iers2010.hpp"
#include "iers2010/cel2ter.hpp"
#include "planetpos.hpp"
#include <cassert>

int main(int argc, char *argv[]) {
  if (argc != 5) {
    fprintf(stderr,
            "Usage: %s [07relativistic_icrf.txt]"
            " [00orbit_icrf.txt] [de421.bsp] [naif0012.tls]\n",
            argv[0]);
    return 1;
  }

  /* Load CSPICE/NAIF Kernels */
  dso::cspice::load_if_unloaded_spk(argv[3]);
  dso::cspice::load_if_unloaded_spk(argv[4]);

  /* parse COST-G (acceleration) results */
  std::vector<costg::CostgAcc> avset;
  if (costg::parse_gravity_field(argv[1], avset)) {
    fprintf(stderr, "ERROR Failed to read input acceleration file %s\n",
            argv[1]);
    return 1;
  }

  /* parse COST-G state, ITRF */
  std::vector<costg::CostgExtState> vstate;
  if (costg::parse_satellite_state(argv[2], vstate)) {
    fprintf(stderr, "ERROR Failed to read input state file %s\n", argv[3]);
    return 1;
  }

  /* loop through reference results ... */
  auto state = vstate.begin();
  auto aset = avset.begin();
  int error = 0;
  while (state != vstate.end()) {
    aset = std::find_if(aset, avset.end(), [&](const costg::CostgAcc &g) {
      return g.gpst == state->gpst;
    });
    if (aset == avset.end()) {
      fprintf(stderr, "ERROR Failed to find matching gravity result!\n");
      return 1;
    }

    /* gps time to TAI */
    const auto tai = costg::gps2tai(state->gpst);

    /* Sun position in ICRF */
    Eigen::Matrix<double, 3, 1> sun;
    if (dso::planet_pos(dso::Planet::SUN, tai.tai2tt(), sun))
      ++error;
    if (error) {
      fprintf(stderr, "ERROR Failed extracting Sun/Moon position!\n");
      return 2;
    }

    /* satellite position w.r.t Earth in GCRS */
    const auto r = state->r();
    /* Position of the Earth with respect to the Sun [m] */
    const auto Res = sun;
    /* Time derivative of the position of the Earth with respect to Sun */
    auto dRes = sun;
    {
      const auto dt = 10; /* seconds */
      auto t2 = tai.tai2tt();
      t2.add_seconds(dt);
      Eigen::Matrix<double, 3, 1> sun2;
      if (dso::planet_pos(dso::Planet::SUN, t2, sun2)) {
        fprintf(stderr, "ERROR Failed extracting Sun/Moon position!\n");
        return 2;
      }
      dRes = (sun2-sun)/dt;
    }
    /* satellite velocity, GCRS */
    const auto v = state->v();

    /* Acceleration from relativity in GCRF  */
    Eigen::Matrix<double, 3, 3> partials;
    Eigen::Matrix<double, 3, 1> set_acc;
    set_acc = iers2010::relativistic_correction(r,Res,dRes,v);

    /* print results */
    printf("%.12e %.20e %.20e %.20e %.20e %.20e %.20e\n", aset->gpst.as_mjd(),
           aset->ax, aset->ay, aset->az, set_acc(0), set_acc(1), set_acc(2));

    /* augment state to next */
    ++state;
  }

  return 0;
}
