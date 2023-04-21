#include "astrodynamics.hpp"
#include "costg_utils.hpp"
#include "egravity.hpp"
#include "orbit_integration.hpp"
#include "planetpos.hpp"
#include "iers2010/cel2ter.hpp"
#include <cassert>

/* extracted from COST-G documentation */
constexpr const double GM_Sun = 1.32712442076e20;
constexpr const double GM_Moon = 0.49028010560e13;

int main(int argc, char *argv[]) {
  if (argc != 7) {
    fprintf(stderr,
            "Usage: %s [04solidEarthTide_icrf.txt] [eopc04_14_IAU2000.62-now]"
            " [00orbit_itrf.txt] [de421.bsp] [naif0012.tls] []\n",
            argv[0]);
    return 1;
  }

  /* Load CSPICE/NAIF Kernels */
  dso::cspice::load_if_unloaded_spk(argv[4]);
  dso::cspice::load_if_unloaded_spk(argv[5]);

  /* parse COST-G (acceleration) results */
  std::vector<costg::CostgAcc> avset;
  if (costg::parse_gravity_field(argv[1], avset)) {
    fprintf(stderr, "ERROR Failed to read input acceleration file %s\n",
            argv[1]);
    return 1;
  }

  /* parse COST-G state, ITRF */
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
  if (costg::parse_rotation_quaternions(argv[6], qvrots)) {
    fprintf(stderr, "ERROR Failed to read input state file %s\n", argv[6]);
    return 1;
  }
#endif

  /* create a SolidEarthTideInstance */
  dso::SolidEarthTide SeTide(iers2010::GMe, iers2010::Re, GM_Moon, GM_Sun);

  /* loop through reference results ... */
  auto state = vstate.begin();
  auto aset = avset.begin();
#ifndef USE_OWN_ROTATION_COSTG
  auto quat = qvrots.begin();
#endif
  int error = 0;
  while (state != vstate.end()) {
    aset= std::find_if(aset, avset.end(), [&](const costg::CostgAcc &g) {
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

    /* Sun and Moon position in ICRF */
    Eigen::Matrix<double, 3, 1> sun,moon;
    if (dso::planet_pos(dso::Planet::SUN, tai.tai2tt(), sun)) ++error;
    if (dso::planet_pos(dso::Planet::MOON, tai.tai2tt(), moon)) ++error;
    if (error) {
      fprintf(stderr, "ERROR Failed extracting Sun/Moon position!\n");
      return 2;
    }

    R.prepare(tai.tai2tt());

    /* rotate, ICRF to ITRF */
#ifdef USE_OWN_ROTATION_COSTG
    sun = R.gcrf2itrf(sun);
    moon = R.gcrf2itrf(moon);
    auto r = R.gcrf2itrf(state->r());
#else
    const Eigen::Quaternion<double> q = quat->q.conjugate();
    sun = q * sun;
    moon = q * moon;
    auto r = q * state->r();
#endif
    

    /* Acceleration from Solid E. T. in ITRF  */
    Eigen::Matrix<double, 3, 3> partials;
    Eigen::Matrix<double, 3, 1> set_acc;
    if (SeTide.acceleration(tai.tai2tt(), R.ut1(), r, moon, sun, set_acc, partials)) {
      fprintf(stderr, "ERROR Failed to compute Solid Earth Tide!\n");
      return 1;
    }

    /* Acceleration from Solid E. T. in GCRF */
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
