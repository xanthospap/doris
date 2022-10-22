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
#include "iers2010/iers2010.hpp"
#include "iers2010/iersc.hpp"
#include "iers2010/tropo.hpp"
#include "integrators.hpp"
#include "satellites/jason3.hpp"
#include "satellites/jason3_quaternions.hpp"
#include "sp3/sp3.hpp"
#include "sp3/sv_interpolate.hpp" /* debug mode */
#include "var_utils.hpp"
#include <cstdio>
#include <cassert>

// usually using these datetimes ...
using Datetime = dso::datetime<dso::nanoseconds>;

struct TropoDetails {
  double Lhz, mfh;
  double Lwz, mfw;
  double sum() const noexcept { return Lhz * mfh + Lwz * mfw; }
  double sum(double Lwz_) const noexcept { return Lhz * mfh + Lwz_ * mfw; }
  void dump() const noexcept {
    printf(": %.6f * %.6f + %.6f * %.6f = %.6f\n", Lhz,
           mfh, Lwz, mfw, sum());
  }
};

int get_tropo(const dso::datetime<dso::nanoseconds> &t,
              const Eigen::Matrix<double, 3, 1> &bxyz, double zd,
              const dso::Gpt3Grid &grid, TropoDetails &Dtrop) noexcept;

int main(int argc, char *argv[]) {
  
  // check input
  if (argc != 2) {
    fprintf(stderr, "USAGE: %s [YAML CONFIG]\n", argv[0]);
    return 1;
  }

  Eigen::Matrix<double,3,1> xyz;
  xyz << 2890643.937260e0, 1310312.136289e0, 5513963.193492e0;

  // resolve the yaml config file and get the root node
  const YAML::Node config = YAML::LoadFile(argv[1]);
  char buf[256];
  
  // Troposphere
  // -------------------------------------------------------------------------
  // read-in the grid file. That is all for now
  dso::get_yaml_value_depth3(config, "troposphere", "gpt3", "grid", buf);
  dso::Gpt3Grid gpt3_grid(buf);
  TropoDetails cDtropo;

  const Datetime t1 = Datetime(dso::year(2021), dso::month(1),
                               dso::day_of_month(1), dso::nanoseconds(0));
  const Datetime t2 = Datetime(dso::year(2022), dso::month(1),
                               dso::day_of_month(1), dso::nanoseconds(0));
  for (auto t = t1; t < t2; t.add_seconds(dso::nanoseconds(30L*1'000'000'000L))) {
    if (get_tropo(t, xyz, dso::deg2rad(30e0), gpt3_grid, cDtropo)) {
      return 3;
    }
    dso::strftime_ymd_hmfs(t, buf);
    printf("%s %.9f ", buf, t.as_mjd());
    cDtropo.dump();
  }

  return 0;
}

int get_tropo(const dso::datetime<dso::nanoseconds> &t,
              const Eigen::Matrix<double, 3, 1> &bxyz, double zd,
              const dso::Gpt3Grid &grid, TropoDetails &Dtrop) noexcept {

  // validate zenith angle
  if (!(zd >= 0e0 && zd <= dso::DPI / 2e0)) {
    fprintf(stderr, "WTF!! Weird zenith angle, is: %.2f\n", dso::rad2deg(zd));
  }
  assert(zd >= 0e0 && zd <= dso::DPI / 2e0);

  // ellipsoidal coordinates of the station; store them in an array
  Eigen::Matrix<double, 3, 1> bell = dso::car2ell<dso::ellipsoid::grs80>(bxyz);
  std::vector<std::array<double, 3>> ellipsoidal(
      1, std::array<double, 3>{bell(0), bell(1), bell(2)});

  // store GPT3 results here
  std::vector<dso::gpt3_result> g3out;

  // call gpt3_fast to get parameters; store them at g3out
  if (dso::gpt3_fast(t, ellipsoidal, 0, grid, g3out)) {
    fprintf(stderr, "[ERROR] Failed to compute gpt3!\n");
    return 20;
  }

  // use VMF3 to compute mapping functions
  double mfh, mfw;
  if (dso::vmf3(g3out[0].ah, g3out[0].aw, t, bell(1), bell(0), zd, mfh, mfw)) {
    fprintf(stderr, "Failed to compute VMF3\n");
    return 30;
  }

  // use refined saastamnoinen to compute the hydrostatic delay (zenith)
  const double zhd0 = dso::saasthyd(g3out[0].p, bell(1), bell(2));
  // apply VMF3 mapping function to compute hydrostatic delay at given zenith
  // const double Dtropo_hydrostatic = zhd0 * vmf3_res.mfh;

  // use Askne and Nordius to approximate wet dealy in zenith
  const double zwd0 = dso::asknewet(g3out[0].e, g3out[0].Tm, g3out[0].la);
  // apply VMF3 mapping function
  // const double Dtropo_wet = zwd0 * vmf3_res.mfw;

  Dtrop.Lhz = zhd0;
  Dtrop.mfh = mfh;
  Dtrop.Lwz = zwd0;
  Dtrop.mfw = mfw;

  return 0;
}
