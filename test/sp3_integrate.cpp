#include "astrodynamics.hpp"
#include "beacon_tbl.hpp"
#include "datetime/datetime_write.hpp"
#include "datetime/utcdates.hpp"
#include "datetime/dtcalendar.hpp"
#include "doris_observation_equations.hpp"
#include "doris_rinex.hpp"
#include "doris_utils.hpp"
#include "filters/filters.hpp"
#include "tides.hpp"
#include "geodesy/geoconst.hpp"
#include "geodesy/geodesy.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iers2010.hpp"
#include "iers2010/iersc.hpp"
#include "iers2010/tropo.hpp"
#include "iers2010/hardisp.hpp"
#include "integrators.hpp"
#include "satellites/jason3.hpp"
#include "satellites/jason3_quaternions.hpp"
#include "sp3/sp3.hpp"
#include "sp3/sv_interpolate.hpp"
#include "var_utils.hpp"
#include <cassert>
#include <cstdio>
#include <cstring>

int main(int argc, char *argv[]) {
  // Check input
  if (argc != 2) {
    fprintf(stderr, "USAGE: %s [YAML CONFIG]\n", argv[0]);
    return 1;
  }

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
    error = dso::parse_gravity_model(buf, degree, order, rnx.ref_datetime(),
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
  dso::IntegrationParameters IntegrationParams(degree, order, eop_lut,
                                               harmonics, buf);

  /*
  IntegrationParams.macromodel =
      dso::MacroModel<dso::SATELLITE::Jason3>::mmcomponents;
  IntegrationParams.numMacroModelComponents =
      dso::MacroModel<dso::SATELLITE::Jason3>::NumPlates;
  IntegrationParams.qhunt = &qhunt;
  IntegrationParams.SatMass = &sat_mass;
  */

  dso::Nrlmsise00 nrlmsise00;
  IntegrationParams.nrlmsise00 = &nrlmsise00;

  dso::get_yaml_value_depth3(config, "force-model", "atmospheric-drag",
                             "atmo-data-csv", buf);
  dso::modified_julian_day utc_mjd_fo;
  const double utc_fday = dso::tai2utc(rnx.time_of_first_obs(), utc_mjd_fo);
  dso::nrlmsise00::InParams<
      dso::nrlmsise00::detail::FluxDataFeedType::ST_CSV_SW>
      atm_data_feed(buf, utc_mjd_fo, utc_fday * 86400e0);
  atm_data_feed.params_.set_switches_on();
  atm_data_feed.params_.use_aparray();
  atm_data_feed.params_.meters_on();
  IntegrationParams.AtmDataFeed = &atm_data_feed;

  // Orbit Integrator
  // -------------------------------------------------------------------------
  // Setup an integrator, to extrapolate orbit with:
  // 1. Relative accuracy 1e-12
  // 2. Absolute accuracy 1e-12
  // 3. Num of Equations: 6 for state and 6*6 for variational equations
  dso::SGOde Integrator(dso::VariationalEquations, 6 + 6 * 6 + 6 * Np, 1e-12,
                        1e-15, &IntegrationParams);
  
  // Important !!
  // set integration parameter estimates to point to the Filter
  // IntegrationParams.estimates = &Filter.x;
  IntegrationParams.drag_coef = &(Filter._ekf._ekf.x(Filter.drag_coef_index()));

  // Important !! OC-TIDES
  IntegrationParams.octide = &octide;

  // Setup Solid Earth Tide
  dso::SolidEarthTide setide(harmonics.GM(), harmonics.Re(),
                             IntegrationParams.GMMon * 1e9,
                             IntegrationParams.GMSun * 1e9);
  IntegrationParams.setide = &setide;

