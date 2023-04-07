#include "astrodynamics.hpp"
#include "beacon_tbl.hpp"
#include "datetime/datetime_write.hpp"
#include "datetime/utcdates.hpp"
#include "doris_observation_equations.hpp"
#include "doris_rinex.hpp"
#include "doris_utils.hpp"
#include "filters/filters.hpp"
#include "geodesy/geoconst.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "iers2010/cel2ter.hpp"
#include "iers2010/hardisp.hpp"
#include "iers2010/iers2010.hpp"
#include "iers2010/iersc.hpp"
#include "iers2010/tropo.hpp"
#include "integrators.hpp"
#include "satellites/jason3.hpp"
#include "satellites/jason3_quaternions.hpp"
#include "sp3/sp3.hpp"
#include "sp3/sv_interpolate.hpp" /* debug mode */
#include "tides.hpp"
#include "var_utils.hpp"
#include <cassert>
#include <cstdio>
#include <cstring>
#include <datetime/dtcalendar.hpp>
#include <datetime/dtfund.hpp>
#include <iers2010/eop.hpp>

constexpr const int INCLUDE_EARTH_TIDES = true;
constexpr const int INCLUDE_OCEAN_TIDES = true;

int integrate(const dso::Sp3DataBlock &sp3block,
              const dso::TwoPartDate &t_target, dso::SGOde &integrator,
              dso::IntegrationParameters &params,
              Eigen::Matrix<double, 6, 1> &state, Eigen::Matrix<double, 6, 1> &state_eci) {
  // reference state from sp3 block (terrestrial to celestial)
  Eigen::Matrix<double, 6, 1> y;
  y << sp3block.state[0] * 1e3, sp3block.state[1] * 1e3,
      sp3block.state[2] * 1e3, sp3block.state[4] * 1e-1,
      sp3block.state[5] * 1e-1, sp3block.state[6] * 1e-1;

  // terrestrial to celestial for reference epoch
  const dso::TwoPartDate t0(sp3block.t);
  dso::Itrs2Gcrs Rot(t0.tai2tt(), &params.eopLUT);
  y = Rot.itrf2gcrf(y);
  const Eigen::Matrix<double,6,1> Y0 = y;

  // variational equations
  Eigen::Matrix<double, 6, 6> I = Eigen::Matrix<double, 6, 6>::Identity();
  I.transposeInPlace();
  Eigen::Matrix<double, 6 * 6, 1> Iv =
      Eigen::Map<Eigen::Matrix<double, 6 * 6, 1>>(I.data(),
                                                  I.cols() * I.rows());
  constexpr const int NumEqns = 6 + (6 * 6);
  static_assert(NumEqns >= 6);
  Eigen::Matrix<double, NumEqns, 1> yF0;
  yF0.block<6 * 6, 1>(6, 0) = Iv;
  yF0.block<6, 1>(0, 0) = y;

  Eigen::VectorXd solution = Eigen::Matrix<double, NumEqns, 1>::Zero();

  double t = 0e0;
  params.mjd_tai = dso::TwoPartDate(t0);
  const auto dt = t_target - params.mjd_tai;
  double tout = dt._big * 86400e0 + dt._small * 86400e0;
  integrator.flag() = dso::SGOde::IFLAG::RESTART;

  fprintf(stderr, "[DEBUG] ODE system of %dx%d equations\n", (int)yF0.rows(),
          (int)yF0.cols());
  if (integrator.de(t, tout, yF0, solution) != dso::SGOde::IFLAG::SUCCESS) {
    fprintf(stderr, "Integration failed!\n");
    return 1;
  }

  state = solution.block<6, 1>(0, 0); // inertial
  dso::TwoPartDate tres(params.mjd_tai);
  tres._small += (tout / 86400e0);
  tres.normalize();
  state_eci = state;

  // check results
  /*{
    const Eigen::Matrix<double,6,1> _y = solution.block<6,1>(0,0);
    const Eigen::Matrix<double,6,6> _F(solution.data()+6);
    const Eigen::Matrix<double,6,1> _y2 = _F * (Y0 -_y);
    printf("[CHK]                               %+.6f %+.6f %+.6f %+.9f %+.9f %+.9f\n",_y(0)-_y2(0),_y(1)-_y2(1),_y(2)-_y2(2),_y(3)-_y2(3),_y(4)-_y2(4),_y(5)-_y2(5));
  }*/

  // celestial to terestrial for epoch
  Rot.prepare(tres.tai2tt());
  state = Rot.gcrf2itrf(state);
  return 0;
}

std::vector<dso::datetime<dso::nanoseconds>> parse_input(const char *fn) {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "Failed to open input file %s\n", fn);
    return std::vector<dso::datetime<dso::nanoseconds>>{};
  }

  std::vector<dso::datetime<dso::nanoseconds>> vec;
  char line[1024];
  while (fin.getline(line, 1024)) {
    try {
      auto t = dso::strptime_ymd_hms<dso::nanoseconds>(line);
      vec.push_back(t);
    } catch (std::exception &) {
      ;
    }
  }

  return vec;
}

int main(int argc, char *argv[]) {
  // Check input
  if (argc != 3) {
    fprintf(stderr, "USAGE: %s [YAML CONFIG] [INPUT FILE]\n", argv[0]);
    return 1;
  }

  auto tvec = parse_input(argv[2]);
  if (!tvec.size())
    return 0;

  // Resolve the yaml config file and get the root node
  const YAML::Node config = YAML::LoadFile(argv[1]);
  char buf[256];
  int error;

  // Load CSPICE/NAIF Kernels (2/3)
  // -------------------------------------------------------------------------
  if (dso::get_yaml_value_depth2(config, "naif-kernels", "spk", buf)) {
    fprintf(stderr, "ERROR failed to find spk kernel\n");
    return 1;
  }
  dso::cspice::load_if_unloaded_spk(buf);
  if (dso::get_yaml_value_depth2(config, "naif-kernels", "lsk", buf)) {
    fprintf(stderr, "ERROR failed to find spk kernel\n");
    return 1;
  }
  dso::cspice::load_if_unloaded_lsk(buf);

  // Reference Orbit
  // -------------------------------------------------------------------------
  dso::get_yaml_value_depth2(config, "data", "sp3", buf);
  dso::Sp3c sp3(buf);
  dso::Sp3Iterator sp3_iterator(sp3);

  // EOP Look Up Table
  // Parse the input EOP data file to create an EopLookUpTable eop_lut
  // -------------------------------------------------------------------------
  dso::EopLookUpTable eop_lut;
  if (dso::get_yaml_value_depth2(config, "eop-info", "eop-file", buf)) {
    fprintf(stderr,
            "ERROR. Failed parsing eop-info/eop-file file from YAML %s\n",
            argv[1]);
    return 1;
  } else {
    const int ref_mjd = sp3.start_epoch().as_mjd();
    const int start = ref_mjd - 10;
    const int end = ref_mjd + 15;
    // parse C04 EOPs
    if (dso::parse_iers_C0420(buf, start, end, eop_lut)) {
      fprintf(stderr, "ERROR. Failed collecting EOP data\n");
      return 1;
    }
    eop_lut.utc2tt();
    eop_lut.regularize(false);
  }

  // Gravity
  // -------------------------------------------------------------------------
  // parse degree, order and the requested gravity model into a
  // HarmonicCoeffs instance. Note that to compute potential we will need
  // Lagrange polynomials (later on)
  int degree, order;
  error = 0;
  error = dso::get_yaml_value_depth2<int>(config, "gravity", "degree", degree);
  error += dso::get_yaml_value_depth2<int>(config, "gravity", "order", order);
  dso::HarmonicCoeffs harmonics(degree);
  if (!error)
    error = dso::get_yaml_value_depth2(config, "gravity", "model", buf);
  // parse the un-normalized harmonic ceofficients from model (gfc format)
  if (!error)
    error = dso::parse_gravity_model(buf, degree, order, sp3.start_epoch(),
                                     harmonics, false);
  if (error) {
    fprintf(stderr, "ERROR Failed handling gravity field model!\n");
    return 1;
  }

  // Setup Integration Parameters for Orbit Integration
  // -------------------------------------------------------------------------
  // We will need the pck (SPICE) kernel for gravitational parameters of Sun
  // and Moon
  if (dso::get_yaml_value_depth2(config, "naif-kernels", "pck", buf)) {
    fprintf(stderr, "ERROR Failed locating NAIF pck kernel\n");
    return 1;
  }
  dso::IntegrationParameters Params(degree, order, eop_lut, harmonics, buf, "data/DTM_2020_F107_Kp.dat");

  // Orbit Integrator
  // -------------------------------------------------------------------------
  // Setup an integrator, to extrapolate orbit with:
  // 1. Relative accuracy 1e-12
  // 2. Absolute accuracy 1e-12
  // 3. Num of Equations: 6 for state and 6*6 for variational equations
  constexpr const int NumEqns = 6 + (6 * 6);
  static_assert(NumEqns >= 6);
  dso::SGOde Integrator(dso::VariationalEquations2, NumEqns, 1e-12, 1e-15,
                        &Params);

  // Setup Solid Earth Tide
  // -------------------------------------------------------------------------
  dso::SolidEarthTide setide(harmonics.GM(), harmonics.Re(), Params.GMMon * 1e9,
                             Params.GMSun * 1e9);
  Params.setide = &setide;
  if (!INCLUDE_EARTH_TIDES)
    Params.setide = nullptr;

  // Ocean Tides Geopotential
  // ------------------------------------------------------------------------
  std::vector<dso::DoodsonOceanTideConstituent> vdds;
  int oc_degree, oc_order;
  {
    dso::get_yaml_value_depth2(config, "ocean-tides", "harmonics", buf);
    error = dso::get_yaml_value_depth2<int>(config, "ocean-tides", "degree",
                                            oc_degree);
    error += dso::get_yaml_value_depth2<int>(config, "ocean-tides", "order",
                                             oc_order);
    if (error || (oc_order > oc_degree)) {
      fprintf(stderr, "Invalid ocean tide degree/order, %d/%d!\n", oc_degree,
              oc_order);
      return 1;
    }
    if (dso::memmap_octide_coefficients(buf, vdds, oc_degree, oc_order, 3,
                                        1e-11)) {
      fprintf(stderr, "Failed reading ocean tide loading file %s!\n", buf);
      return 1;
    }
  }
  // An ocean tide instance
  dso::OceanTide octide(vdds, harmonics.GM(), harmonics.Re(), oc_degree,
                        oc_order);
  Params.octide = &octide;
  if (!INCLUDE_OCEAN_TIDES)
    Params.octide = nullptr;

  sp3_iterator.begin();
  dso::Sp3DataBlock sp3block;
  auto tend = std::unique(tvec.begin(), tvec.end());
  error = 0;
  for (auto t = tvec.begin(); t != tend; ++t) {
    dso::strftime_ymd_hmfs<dso::nanoseconds>(*t, buf);
    if (sp3_iterator.goto_epoch(*t)) {
      fprintf(stderr, "ERROR Failed to get reference position from SP3\n");
      return 1;
    }
    // assign current block
    sp3block = sp3_iterator.data_block();
    // reference time
    const dso::TwoPartDate t0(sp3block.t);
    // integrate
    Eigen::Matrix<double, 6, 1> state, state_eci;
    if (integrate(sp3block, dso::TwoPartDate(*t), Integrator, Params, state, state_eci)) {
      return 2;
    }
    // print results in terrestrial RF
    dso::strftime_ymd_hmfs<dso::nanoseconds>(*t, buf);
    printf("%s %+.6f %+.6f %+.6f %+.9f %+.9f %+.9f\n", buf, state(0),
           state(1), state(2), state(3), state(4), state(5));
    printf("[ECI] %s %+.9f %+.9f %+.9f %+.12e %+.12e %+.12e\n", buf, state_eci(0),
           state_eci(1), state_eci(2), state_eci(3), state_eci(4),
           state_eci(5));
    const dso::TwoPartDate t1(*t);
    printf("#[CMT] Integration interval: %.3f[sec]\n",
           t1.diff<dso::DateTimeDifferenceType::FractionalSeconds>(t0));
  }

  return 0;
}
