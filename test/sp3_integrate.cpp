#include "astrodynamics.hpp"
#include "beacon_tbl.hpp"
#include "datetime/datetime_write.hpp"
#include "datetime/dtcalendar.hpp"
#include "datetime/utcdates.hpp"
#include "doris_observation_equations.hpp"
#include "doris_rinex.hpp"
#include "doris_utils.hpp"
#include "filters/filters.hpp"
#include "geodesy/geoconst.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "iers2010/hardisp.hpp"
#include "iers2010/iers2010.hpp"
#include "iers2010/iersc.hpp"
#include "iers2010/tropo.hpp"
#include "integrators.hpp"
#include "satellites/jason3.hpp"
#include "satellites/jason3_quaternions.hpp"
#include "sp3/sp3.hpp"
#include "sp3/sv_interpolate.hpp"
#include "tides.hpp"
#include "var_utils.hpp"
#include <cassert>
#include <cstdio>
#include <cstring>
#include <datetime/dtfund.hpp>

constexpr const double MAXHOURS = 48e0;
constexpr const int INCLUDE_EARTH_TIDES = true;

dso::datetime<dso::nanoseconds> dttr(const dso::TwoPartDate &t) {
  const auto t1 = t.normalized();
  const long inano = static_cast<long>((t1._small * 86400e0) *
                                       dso::nanoseconds::sec_factor<double>());
  return dso::datetime<dso::nanoseconds>(t1._big, dso::nanoseconds(inano));
}

int integrate(const dso::Sp3DataBlock &sp3block,
              const dso::TwoPartDate &t_target, dso::SGOde &integrator,
              dso::IntegrationParameters &params,
              Eigen::Matrix<double, 6, 1> &state) {
  // reference state from sp3 block (terrestrial to celestial)
  Eigen::Matrix<double, 6, 1> y;
  y << sp3block.state[0] * 1e3, sp3block.state[1] * 1e3,
      sp3block.state[2] * 1e3, sp3block.state[4] * 1e-1,
      sp3block.state[5] * 1e-1, sp3block.state[6] * 1e-1;

  // terrestrial to celestial for reference epoch
  Eigen::Matrix<double, 3, 3> rc2i, rpom;
  double era, xlod;
  assert(!gcrs2itrs(dso::TwoPartDate(sp3block.t), params.eopLUT, rc2i, era,
                    rpom, xlod));
  y = dso::yter2cel(y, rc2i, era, xlod, rpom);

  // variational equations
  Eigen::Matrix<double, 6, 6> I = Eigen::Matrix<double, 6, 6>::Identity();
  I.transposeInPlace();
  Eigen::Matrix<double, 6 * 6, 1> Iv =
      Eigen::Map<Eigen::Matrix<double, 6 * 6, 1>>(I.data(),
                                                  I.cols() * I.rows());
  Eigen::Matrix<double, 6 + 6 * 6, 1> yF0;
  yF0.block<6, 1>(0, 0) = y;
  yF0.block<6 * 6, 1>(6, 0) = Iv;

  // Eigen::Matrix<double, 6 + 6 * 6, 1> solution;
  Eigen::VectorXd solution = Eigen::Matrix<double,6+6*6,1>::Zero();

  double t = 0e0;
  params.mjd_tai = dso::TwoPartDate(sp3block.t);
  const auto dt = t_target - params.mjd_tai;
  double tout = dt._big * 86400e0 + dt._small * 86400e0;
  integrator.flag() = dso::SGOde::IFLAG::RESTART;

  if (integrator.de(t, tout, yF0, solution) != dso::SGOde::IFLAG::SUCCESS) {
    fprintf(stderr, "Integration failed!\n");
    return 1;
  }

  state = solution.block<6, 1>(0, 0); // inertial
  dso::TwoPartDate tres(params.mjd_tai);
  tres._small += (tout / 86400e0);
  tres.normalize();
  // terrestrial to celestial for reference epoch
  assert(!gcrs2itrs(tres, params.eopLUT, rc2i, era, rpom, xlod));
  state = dso::ycel2ter(state, rc2i, era, xlod, rpom);
  return 0;
}

int main(int argc, char *argv[]) {
  // Check input
  if (argc != 2 && argc!= 3) {
    fprintf(stderr, "USAGE: %s [YAML CONFIG] <INTEGRATION_MIN>\n", argv[0]);
    return 1;
  }

  int INTEGRATION_MIN = 30;
  if (argc == 3) {
    INTEGRATION_MIN = std::atoi(argv[2]);
    assert(INTEGRATION_MIN>0);
  }
  double IDAYS = (INTEGRATION_MIN / 24e0 / 60e0);

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
    const int start = ref_mjd - 5;
    const int end = ref_mjd + 10;
    // parse C04 EOPs
    if (parse_iers_C04(buf, start, end, eop_lut)) {
      fprintf(stderr, "ERROR. Failed collecting EOP data\n");
      return 1;
    }
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
  dso::IntegrationParameters Params(degree, order, eop_lut, harmonics, buf);

  // Orbit Integrator
  // -------------------------------------------------------------------------
  // Setup an integrator, to extrapolate orbit with:
  // 1. Relative accuracy 1e-12
  // 2. Absolute accuracy 1e-12
  // 3. Num of Equations: 6 for state and 6*6 for variational equations
  dso::SGOde Integrator(dso::VariationalEquations, 6 + 6 * 6, 1e-12, 1e-15,
                        &Params);

  // Setup Solid Earth Tide
  // -------------------------------------------------------------------------
  dso::SolidEarthTide setide(harmonics.GM(), harmonics.Re(),
                             Params.GMMon * 1e9,
                             Params.GMSun * 1e9);
  Params.setide = &setide;
  if (!INCLUDE_EARTH_TIDES) Params.setide = nullptr;

  sp3_iterator.begin();
  const dso::TwoPartDate t_start =
      dso::TwoPartDate(sp3_iterator.data_block().t);
  dso::Sp3DataBlock sp3block;
  error = 0;
  long count_epochs = 0;
  while (!error) {
    sp3block = sp3_iterator.data_block();
    // reference time
    const dso::TwoPartDate t0(sp3block.t);
    // print results in terrestrial RF
    dso::strftime_ymd_hmfs(dttr(t0), buf);
    printf("[SP3] %s %.12e %.12e %.12e %.12e %.12e %.12e\n", buf,
           sp3block.state[0] * 1e3, sp3block.state[1] * 1e3,
           sp3block.state[2] * 1e3, sp3block.state[4] * 1e-1,
           sp3block.state[5] * 1e-1, sp3block.state[6] * 1e-1);
    // target integration time
    const dso::TwoPartDate tt(t0._big, t0._small + IDAYS);
    const dso::TwoPartDate t = tt.normalized();
    // integrate
    if (!(count_epochs % 5)) {
      Eigen::Matrix<double, 6, 1> state;
      if (integrate(sp3block, t, Integrator, Params, state)) {
        return 2;
      }
      // print results in terrestrial RF
      dso::strftime_ymd_hmfs(dttr(t), buf);
      printf("[INT] %s %.12e %.12e %.12e %.12e %.12e %.12e\n", buf,
             state(0), state(1), state(2), state(3), state(4), state(5));
    }
    if (t0-t_start > dso::TwoPartDate(0,MAXHOURS/24e0)) break;
    // next sp3 epoch ...
    if (sp3_iterator.advance()) {
      fprintf(stderr, "ERROR Failed reading next epoch from sp3 file!\n");
      return 10;
    }
    ++count_epochs;
  }

  return 0;
}
