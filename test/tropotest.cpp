#include "iers2010/tropo.hpp"
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
#include "integrators.hpp"
#include "satellites/jason3.hpp"
#include "satellites/jason3_quaternions.hpp"
#include "sp3/sp3.hpp"
#include "sp3/sv_interpolate.hpp" /* debug mode */
#include "var_utils.hpp"
#include <cassert>
#include <cstdio>

// usually using these datetimes ...
using Datetime = dso::datetime<dso::nanoseconds>;

struct TropoDetails {
  double Lhz, mfh;
  double Lwz, mfw;
  double sum() const noexcept { return Lhz * mfh + Lwz * mfw; }
  double sum(double Lwz_) const noexcept { return Lhz * mfh + Lwz_ * mfw; }
  void dump(FILE *fp, double el_deg) const noexcept {
    fprintf(fp, " %.6f %.6f %.6f %.6f %.6f %.1f\n", Lhz, mfh, Lwz, mfw, sum(),
            el_deg);
  }
  void dump(double el_deg) const noexcept {
    printf(" %.6f %.6f %.6f %.6f %.6f %.1f\n", Lhz, mfh, Lwz, mfw, sum(),
           el_deg);
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

  // Eigen::Matrix<double,3,1> xyz;
  // xyz << 2890643.937260e0, 1310312.136289e0, 5513963.193492e0;
  std::vector<Eigen::Matrix<double, 3, 1>> xyzs(6);
  // TERRE-ADELIE - ANTARCTICA,
  // https://ids-doris.org/network/sitelogs/station.html?code=TERRE-ADELIE
  xyzs[0] << -1940878.515e0, 1628473.041e0, -5833723.413e0;
  //// NY-ALESUND II - NORWAY,
  /// https://ids-doris.org/network/sitelogs/station.html?code=NY-ALESUND%20II
  /// xyzs[1] << 1201300.046e0, 251874.435e0, 6238000.324e0;
  
  /// DIONYSOS - GREECE,
  /// https://ids-doris.org/network/sitelogs/station.html?code=DIONYSOS
  xyzs[1] << 4595212.468e0, 2039473.691e0, 3912617.891e0;
  
  /// https://apps.ids-doris.org/server/php/api/stations/doris/HBMB/sitelog/HBMB20171024.LOG
  xyzs[2] << 5084653.440e0, 2670347.412e0, -2768470.902e0;
  
  /// KOUROU - FRANCE (French Guiana),
  /// https://ids-doris.org/network/sitelogs/station.html?code=KOUROU
  xyzs[3] << 3855260.440e0, -5049735.535e0, 563056.590e0;
  
  /// https://apps.ids-doris.org/server/php/api/stations/doris/THUB/sitelog/THUB20170301.LOG
  xyzs[4] << 538110.515e0, -1389031.364e0, 6180994.514e0;
  
  /// https://ids-doris.org/network/sitelogs/station.html?code=TOULOUSE
  xyzs[5] << 4628693.773e0, 119984.893e0, 4372104.604e0;
  /// ASCENSION - UNITED KINGDOM,
  /// https://ids-doris.org/network/sitelogs/station.html?code=ASCENSION
  /// xyzs[4] << 6121154.081e0, -1563976.723e0, -872606.019e0;
  /// WETTZELL - GERMANY,
  /// https://ids-doris.org/network/sitelogs/station.html?code=WETTZELL
  /// xyzs[5] << 4075559.776e0, 931580.241e0, 4801621.584e0;
  [[maybe_unused]] const char *names[] = {"ADHC", "DIOB", "HBMB", "KRWB", "THUB", "TLSB"};

  // resolve the yaml config file and get the root node
  const YAML::Node config = YAML::LoadFile(argv[1]);
  char buf[256];

  // Troposphere
  // -------------------------------------------------------------------------
  // read-in the grid file. That is all for now
  dso::get_yaml_value_depth3(config, "troposphere", "gpt3", "grid", buf);
  dso::Gpt3Grid gpt3_grid(buf);
  [[maybe_unused]] TropoDetails cDtropo;

  int it = 0;
  for (const Eigen::Matrix<double, 3, 1> &xyz : xyzs) {

    // ellipsoidal coordinates of the station; store them in an array
    Eigen::Matrix<double, 3, 1> ell =
        dso::car2ell<dso::ellipsoid::grs80>(xyz);
    std::vector<std::array<double, 3>> ellipsoidal(
        1, std::array<double, 3>{ell(0), ell(1), ell(2)});

    // for (int elevation = 10; elevation < 90; elevation += 10) {
      
      // open output file
      FILE *fp;
      char fout[64];
      std::strcpy(fout, names[it]);
      std::strcat(fout, "-gpt.dat");
      // sprintf(fout + std::strlen(fout), "%02d", elevation);

      fp = fopen(fout, "w");
      printf("Writing output to file: %s\n", fout);
      //fprintf(fp, "%s %s %s %s %s %s %s %s %s\n", "Date", "Time", "Mjd", "Lhz", "mfh",
      //        "Lwz", "mfw", "Dtrop", "Elevation");
      // fprintf(fp, "%s %s %s %s %s %s %s\n", "Mjd", "Lhz", "mfh",
      //         "Lwz", "mfw", "Dtrop", "Elevation");

      const Datetime t1 = Datetime(dso::year(2021), dso::month(1),
                                   dso::day_of_month(1), dso::nanoseconds(0));
      const Datetime t2 = Datetime(dso::year(2022), dso::month(1),
                                   dso::day_of_month(1), dso::nanoseconds(0));
      for (auto t = t1; t < t2;
           t.add_seconds(dso::nanoseconds(7L * 1'000'000'000L))) {
        // const double zd = dso::deg2rad(90 - elevation);
        // if (get_tropo(t, xyz, zd, gpt3_grid, cDtropo)) {
        //   return 3;
        // }

        // store GPT3 results here
        std::vector<dso::gpt3_result> g3out;

        // call gpt3_fast to get parameters; store them at g3out
        if (dso::gpt3_fast(t, ellipsoidal, 0, gpt3_grid, g3out)) {
          fprintf(stderr, "[ERROR] Failed to compute gpt3!\n");
          return 20;
        }

        dso::strftime_ymd_hmfs(t, buf);
        const double es =  6.11e0 * std::exp(17.27e0 * g3out[0].T / (237.3e0 + g3out[0].T));
        fprintf(fp, "%s, %s, %.3f, %.3f, %.3f, %.3f, %.9f\n", buf, names[it],
                g3out[0].p, g3out[0].T, (g3out[0].e / es) * 100e0,
                g3out[0].e, t.as_mjd());
        // cDtropo.dump(fp, (double)elevation);
      }

      fclose(fp);
    // }
    ++it;
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
