#include "costg_utils.hpp"
#include "astrodynamics.hpp"
#include "egravity.hpp"
#include <cassert>

constexpr const int degree = 180;
constexpr const int order = 180;

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr,
            "Usage: %s [GRAVITY FIELD] [02gravityfield_itrf.txt] "
            "[00orbit_itrf.txt]\n",
            argv[0]);
    return 1;
  }

  /* parse COST-G (acceleration) results */
  std::vector<costg::CostgAcc> vgrav;
  if (costg::parse_gravity_field(argv[2], vgrav)) {
    fprintf(stderr, "ERROR Failed to read input acceleration file %s\n",
            argv[2]);
    return 1;
  }
  
  /* parse COST-G state */
  std::vector<costg::CostgExtState> vstate;
  if (costg::parse_satellite_state(argv[3], vstate)) {
    fprintf(stderr, "ERROR Failed to read input state file %s\n",
            argv[3]);
    return 1;
  }

  dso::StokesCoeffs stokes(degree);
  if (dso::parse_gravity_model(argv[1], degree, order, vgrav[0].gpst, stokes)) {
    fprintf(stderr, "ERROR Failed parsing gravity model %s\n", argv[1]);
    return 1;
  }

  /* loop through reference results ... */
  auto state = vstate.begin();
  auto rgrav = vgrav.begin();
  while (state != vstate.end()) {
    rgrav = std::find_if(rgrav, vgrav.end(), [&](const costg::CostgAcc& g){return g.gpst == state->gpst;});
    if (rgrav == vgrav.end()) {
      fprintf(stderr, "ERROR Failed to find matching gravity result!\n");
      return 1;
    }

    /* gps time to TAI */
    // const auto tai = costg::gps2tai(rgrav->gpst);

    /* my result */
    Eigen::Matrix<double, 3, 3> gpartials;
    Eigen::Matrix<double, 3, 1> gacc, sphericalacc;
    
    dso::gravity_acceleration(stokes, state->r(), degree, stokes.Re(),
                              stokes.GM(), gacc, gpartials);
    /* compute spherical term degree = 0 */
    dso::gravity_acceleration(stokes, state->r(), 1, stokes.Re(), stokes.GM(),
                              sphericalacc, gpartials);
    gacc -= sphericalacc;

    /* print results */
    printf("%.12e %.20e %.20e %.20e %.20e %.20e %.20e\n", rgrav->gpst.as_mjd(),
           rgrav->ax, rgrav->ay, rgrav->az, gacc(0), gacc(1), gacc(2));

    /* augment state to next */
    ++state;
  }


  return 0;
}
