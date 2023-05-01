#include "astrodynamics.hpp"
#include "costg_utils.hpp"
#include "egravity.hpp"
#include "iers2010/cel2ter.hpp"
#include "orbit_integration.hpp"
#include "planetpos.hpp"
#include <cassert>

constexpr const int degree = 180;
constexpr const int order = 180;

int main(int argc, char *argv[]) {
  if (argc != 6) {
    fprintf(stderr,
            "Usage: %s [11oceanTide_fes2014b_34major_icrf.txt] [eopc04_14_IAU2000.62-now]"
            " [00orbit_icrf.txt] [oceanTide_FES2014b.potential.iers.txt] "
            "[01earthRotation_quaternion.txt]\n",
            argv[0]);
    return 1;
  }

  /* parse COST-G (acceleration) results */
  std::vector<costg::CostgAcc> avset;
  if (costg::parse_gravity_field(argv[1], avset)) {
    fprintf(stderr, "ERROR Failed to read input acceleration file %s\n",
            argv[1]);
    return 1;
  }

  /* parse COST-G state, ICRF */
  std::vector<costg::CostgExtState> vstate;
  if (costg::parse_satellite_state(argv[3], vstate)) {
    fprintf(stderr, "ERROR Failed to read input state file %s\n", argv[3]);
    return 1;
  }
  
  /* create and fill an EOP look up table */
  dso::EopLookUpTable eop_lut;
  {
    const int ref_mjd = vstate[0].gpst.as_mjd();
    const int start = ref_mjd - 3;
    const int end = ref_mjd + 4;
    /* parse C04 EOPs */
    if (parse_iers_C0414(argv[2], start, end, eop_lut)) {
      fprintf(stderr, "ERROR. Failed collecting EOP data\n");
      return 1;
    }
    /* transform the tabulated values (time) to the TT timescale */
    eop_lut.utc2tt();
    /* regularize ΔUT1 and LOD values (to ΔUT1-R and LOD-R) */
    eop_lut.regularize(false);
  }

  /* an ITRF-to-GCRF Rotation for further use */
  dso::Itrs2Gcrs R(costg::gps2tai(vstate[0].gpst).tai2tt(), &eop_lut);

#ifndef USE_OWN_ROTATION_COSTG
  /* read rotation quaternions from costg file */
  std::vector<costg::CostgQuat> qvrots;
  if (costg::parse_rotation_quaternions(argv[5], qvrots)) {
    fprintf(stderr, "ERROR Failed to read input state file %s\n", argv[6]);
    return 1;
  }
#endif

  /* create an Ocean Tide instance */
  dso::OceanTide OcTide(argv[4], 1e-11, degree, order, iers2010::GMe,
                        iers2010::Re);

  /* loop through reference results ... */
  auto state = vstate.begin();
  auto aset = avset.begin();
#ifndef USE_OWN_ROTATION_COSTG
  auto quat = qvrots.begin();
#endif
  while (state != vstate.end()) {
    aset = std::find_if(aset, avset.end(), [&](const costg::CostgAcc &g) {
      return g.gpst == state->gpst;
    });
    if (aset == avset.end()) {
      fprintf(stderr, "ERROR Failed to find matching gravity result!\n");
      return 1;
    }
#ifndef USE_OWN_ROTATION_COSTG
    quat = std::find_if(quat, qvrots.end(), [&](const costg::CostgQuat &g) {
      return g.gpst == state->gpst;
    });
    if (quat == qvrots.end()) {
      fprintf(stderr, "ERROR Failed to find matching quaternion result!\n");
      return 1;
    }
#endif

    /* gps time to TAI */
    const auto tai = costg::gps2tai(state->gpst);

    R.prepare(tai.tai2tt());

    /* rotate, ICRF to ITRF */
#ifdef USE_OWN_ROTATION_COSTG
    auto r = R.gcrf2itrf(state->r());
#else
    const Eigen::Quaternion<double> q = quat->q.conjugate();
    auto r = q * state->r();
#endif

    /* Acceleration from Ocean Tides in ITRF  */
    Eigen::Matrix<double, 3, 3> partials;
    Eigen::Matrix<double, 3, 1> set_acc;
    if (OcTide.acceleration(tai.tai2tt(), r, set_acc, partials)) {
      fprintf(stderr, "ERROR Failed to compute Ocean Tide!\n");
      return 1;
    }

    {
      /* compute and subtract degrees [0,1] */
      Eigen::Matrix<double, 3, 1> set01_acc;
      if (OcTide.acceleration(tai.tai2tt(), r, set01_acc, partials, 1, 1)) {
        fprintf(stderr, "ERROR Failed to compute Ocean Tide!\n");
        return 1;
      }
      set_acc = set_acc - set01_acc;
    }

    /* Acceleration in GCRF */
#ifdef USE_OWN_ROTATION_COSTG
    set_acc = R.itrf2gcrf(set_acc);
#else
    set_acc = q.conjugate() * set_acc;
#endif

    /* print results */
    printf("%.12e %.20e %.20e %.20e %.20e %.20e %.20e\n", aset->gpst.as_mjd(),
           aset->ax, aset->ay, aset->az, set_acc(0), set_acc(1), set_acc(2));

    /* augment state to next */
    ++state;
  }

  return 0;
}
