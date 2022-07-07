#include "astrodynamics.hpp"
#include "datetime/datetime_write.hpp"
#include "egravity.hpp"
#include "eigen3/Eigen/Eigen"
#include "geodesy/units.hpp"
#include "icgemio.hpp"
#include "iers2010/iau.hpp"
#include "iers2010/iers2010.hpp"
#include "iers_bulletin.hpp"
#include "integrators.hpp"
#include "sp3/sp3.hpp"
#include <cstdio>
#include <datetime/dtfund.hpp>

using dso::sp3::SatelliteId;
constexpr const int Degree = 20;
constexpr const int Order = 20;
const int Integrate = false;

// approximate number of data points in a bulletin B file (disregarding
// preliminery values)
constexpr const int BBSZ = 40;

using Datetime = dso::datetime<dso::nanoseconds>;

// xp, yp in milliarcsecond [mas] dut1 in millisecond [ms]
int getEop(const Datetime &t, double &xp, double &yp, double &dut1,
           const double *mjd, const double *xpa, const double *ypa,
           const double *ut1a, int bbblocks_size) {
  // apply corrections; first use INTERP to interpolate x,y,ut1 values for
  // the given date (in [mas], [msec])
  if (iers2010::interp_pole(mjd, xpa, ypa, ut1a, bbblocks_size, t.as_mjd(), xp,
                            yp, dut1)) {
    fprintf(stderr, "ERROR. Failed call to interp_pole\n");
    return 1;
  }
  // account for variations in polar motion (Dx,Dy) ocean-tides; results in
  // [μas]
  double dxp, dyp, dut2;
  if (iers2010::ortho_eop(t, dxp, dyp, dut2)) {
    fprintf(stderr, "ERROR. Failed call to ortho_eop!\n");
    return 4;
  }
  xp += dxp * 1e-3;
  yp += dyp * 1e-3;
  dut1 += dut2 * 1e-3;
  // account for libration effects; results in [μas]
  if (iers2010::pmsdnut2(t, dxp, dyp)) {
    fprintf(stderr, "ERROR Failed call to pmsdnut2!\n");
    return 4;
  }
  xp += dxp * 1e-3;
  yp += dyp * 1e-3;

  return 0;
}

int initEopLookup(const char *bulletin_fn, double *mjd, double *xpa,
                  double *ypa, double *ut1a) {
  // parse the Bulletin B file (diregard preliminary values, only use final)
  dso::IersBulletinB_Section1Block bbblocks[35];
  dso::IersBulletinB bulb(bulletin_fn);
  int bbblocks_size = bulb.parse_section1(bbblocks, false);
  if (bbblocks_size <= 0) {
    fprintf(stderr, "Failed parsing Bulletin B file %s\n", bulletin_fn);
    return 2;
  }

  for (int i = 0; i < bbblocks_size; i++)
    mjd[i] = static_cast<double>(bbblocks[i].mjd);
  for (int i = 0; i < bbblocks_size; i++)
    xpa[i] = bbblocks[i].x;
  for (int i = 0; i < bbblocks_size; i++)
    ypa[i] = bbblocks[i].y;
  for (int i = 0; i < bbblocks_size; i++)
    ut1a[i] = bbblocks[i].dut1;

  return bbblocks_size;
}

int getC04(const Datetime &t, char *bulletin_fn) {
  // get the corresponding bulletin b file from IERS (download)
  int error = dso::download_iers_bulletinb_for(t.mjd().as_underlying_type(),
                                               bulletin_fn);
  if (error) {
    fprintf(stderr, "Failed to download Bulletin B file for date: %ld\n",
            t.mjd().as_underlying_type());
    return 1;
  }

  return 0;
}

// to transfer parameters for Variational Equations
struct AuxParams {
  dso::HarmonicCoeffs *hc;
  dso::Mat2D<dso::MatrixStorageType::Trapezoid> *V, *W;
  const Eigen::Matrix<double, 3, 3> *cel2ter;
  int degree, order;
};

// handle gravity filed
int gravity(const char *gfn, int degree, int order, dso::HarmonicCoeffs &hc) {
  dso::Icgem gfc(gfn);

  // parse the header ...
  if (gfc.parse_header()) {
    fprintf(stderr, "ERROR! Failed to parse icgem header!\n");
    return 1;
  }

  assert(degree <= gfc.degree() && order <= degree && hc.degree() == degree);

  // parse data; store coefficients to hc
  if (gfc.parse_data(degree, order, &hc)) {
    fprintf(stderr, "ERROR! Failed to parse harmonic coefficients\n");
    return 1;
  }

  // de-normalize the harmonics coeffs
  hc.denormalize();

  return 0;
}

Eigen::Matrix<double, 3, 3> ter2cel(const dso::datetime<dso::nanoseconds> &tai,
                                    double xp, double yp,
                                    double dut1) noexcept {
  auto tt = tai;
  tt.add_seconds(dso::nanoseconds(32'184'000'000));
  const double mjd_tt = tt.as_mjd();

  dso::nanoseconds leap_seconds(dso::dat(mjd_tt) * 1'000'000'000L);
  auto utc = tai;
  utc.remove_seconds(leap_seconds);

  auto ut1 = utc;
  dso::nanoseconds dut1_ns(
      static_cast<long>(dut1 * 1e6)); // milliseconds to nanoseconds
  ut1.add_seconds(dut1_ns);
  double mjd_ut1 = ut1.as_mjd();

  auto mat =
      iers2010::sofa::c2t06a(dso::mjd0_jd, mjd_tt, dso::mjd0_jd, mjd_ut1,
                             iers2010::DMAS2R * xp, iers2010::DMAS2R * yp);
  // return Eigen::Matrix<double, 3, 3>(mat.data);
  Eigen::Matrix<double, 3, 3> t2c;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      t2c(j, i) = mat(i, j);
  
  return t2c;
}

void VariationalEquations(
    [[maybe_unused]] double t,
    // state and state transition matrix
    const Eigen::VectorXd &yPhi,
    // state derivative and state transition matrix derivative
    Eigen::Ref<Eigen::VectorXd> yPhiP,
    //
    void *pAux) noexcept {
  static unsigned call_nr = 0;

  // void pointer to AuxParameters
  AuxParams *params = static_cast<AuxParams *>(pAux);

  // split position and velocity vectors
  Eigen::Matrix<double, 3, 1> r = yPhi.block<3, 1>(0, 0);
  Eigen::Matrix<double, 3, 1> v = yPhi.block<3, 1>(3, 0);

  // compute gravity-induced acceleration and gradient
  Eigen::Matrix<double, 3, 1> r_geo = *(params->cel2ter) * r;
  Eigen::Matrix<double, 3, 3> gpartials;
  Eigen::Matrix<double, 3, 1> gacc = dso::grav_potential_accel(
      r_geo, params->degree, params->order, *(params->V), *(params->W),
      *(params->hc), gpartials);

  // State transition (skip first column which is the state vector)
  Eigen::Matrix<double, 6, 6> Phi(yPhi.data() + 6);

  // derivative of state transition matrix, aka
  // | v (3x3)   0 (3x3)     I (3x3)   |
  // | a (3x1) da/dr (3x3) da/dv (3x3) |
  // because dv/dr = 0 and
  //         da/dr = I
  Eigen::Matrix<double, 6, 6> dfdy;
  dfdy.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Zero();
  dfdy.block<3, 3>(0, 3) = Eigen::Matrix<double, 3, 3>::Identity();
  dfdy.block<3, 3>(3, 0) = gpartials;
  dfdy.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Zero();

  // Derivative of combined state vector and state transition matrix
  // dPhi = dfdy * Phi;
  Eigen::Matrix<double, 6, 7> yPhip;
  yPhip.block<6, 6>(0, 1) = dfdy * Phi;

  // state derivative (aka [v,a]), in one (first) column
  yPhip.block<3, 1>(0, 0) = v;
  yPhip.block<3, 1>(3, 0) = gacc;

  // matrix to vector (column-wise)
  yPhiP = Eigen::VectorXd(
      Eigen::Map<Eigen::VectorXd>(yPhip.data(), yPhip.cols() * yPhip.rows()));

  ++call_nr;
  return;
}

int main(int argc, char *argv[]) {

  // handle gravity field
  if (argc < 3 || argc > 5) {
    fprintf(stderr,
            "Usage: %s <SP3 FILE> <GRAVITY MODEL FILE> [DEGREE - optional] "
            "[ORDER - optional]\n",
            argv[0]);
    return 1;
  }

  // parse degree and order of gravity field
  int degree = Degree;
  int order = Order;
  if (argc >= 4)
    degree = std::atoi(argv[3]);
  if (argc == 5)
    order = std::atoi(argv[4]);

  // Harmonic coefficients
  printf("* setting up harmonic coefficients ...\n");
  dso::HarmonicCoeffs hc(degree);
  if (gravity(argv[2], degree, order, hc)) {
    fprintf(stderr,
            "ERROR. Failed to compute harmonic coefficients. Aborting ...\n");
    return 1;
  }

  // Lagrange polynomials (depend on position) N+1 (for zero offset) and +2
  // because we are computing potential partials. Allocate space ...
  dso::Mat2D<dso::MatrixStorageType::Trapezoid> V(degree + 3, order + 3),
      W(degree + 3, order + 3);

  // Handle the Sp3 file, prepare for parsing ...
  printf("* handling input sp3 file ...\n");
  dso::Sp3c sp3(argv[1]);
  //#ifdef DEBUG
  // sp3.print_members();
  //#endif

  SatelliteId sv;
  dso::Sp3DataBlock block;

  if (sp3.num_sats() == 1) {
    sv.set_id(sp3.sattellite_vector()[0].id);
  } else {
    fprintf(stderr, "More than one Satellites in sp3 file. This test program "
                    "is only meant to work with one.\n");
    return 1;
  }

  printf("* downloading and parsing Bulletin B file ...\n");
  const auto t0 = sp3.start_epoch();
  char bulletin_fn[256];
  if (getC04(t0, bulletin_fn)) {
    return 3;
  }
  printf("\t-> downloaded...\n");
  double mjd[BBSZ], xpa[BBSZ], ypa[BBSZ], ut1a[BBSZ];
  int bbblocks_size = initEopLookup(bulletin_fn, mjd, xpa, ypa, ut1a);
  printf("\t-> lookup table set up...\n");
  double xp, yp, dut1;
  if (getEop(t0, xp, yp, dut1, mjd, xpa, ypa, ut1a, bbblocks_size)) {
    return 6;
  }
  printf("\t-> got eop data for date!\n");
  printf("\txp=%.15f yp=%.15f dut1=%.15f [mas], [msec]\n", xp, yp, dut1);

  // Create Auxiliary parameters
  AuxParams params{&hc, &V, &W, nullptr, degree, order};
  // SetUp an integrator
  printf("* setting up integrator ...\n");
  const double relerr = 1.0e-9;  // Relative and absolute
  const double abserr = 1.0e-16; // accuracy requirement
  dso::SGOde Integrator(VariationalEquations, 6 * 7, relerr, abserr, &params);
  Eigen::Matrix<double, 6 * 7, 1> yPhi;
  Eigen::VectorXd sol(6 * 7);
  Eigen::VectorXd r0_geo(6), r0_cel(6);

  // Time
  // const auto t0 = sp3.start_epoch();
  auto epoch = t0;                  // current
  const auto step = sp3.interval(); // dso::nanoseconds(1'000'000'000L * 5*60);

  // let's try reading the records; note that -1 denotes EOF
  int error = 0;
  std::size_t rec_count = 0;
  char buf[64];
  printf("* iterating ...\n");
  do {

    // read nex record ...
    if ((error = sp3.get_next_data_block(sv, block)) > 0) {
      printf("Something went wrong ....status = %3d\n", error);
      return 1;
    }

    // check the heath status
    bool position_ok = !block.flag.is_set(dso::Sp3Event::bad_abscent_position);

    if (position_ok) {

      // transform geocentric state to inertial
      Eigen::Matrix<double, 3, 3> t2c(ter2cel(block.t, xp, yp, dut1));
      Eigen::Matrix<double, 3, 3> c2t = t2c.transpose();
      params.cel2ter = &c2t;

      // Vector containing state + variational equations
      yPhi = Eigen::Matrix<double, 6 * 7, 1>::Zero();

      // accumulate state (m, m/sec)
      r0_geo << block.state[0] * 1e3, block.state[1] * 1e3,
          block.state[2] * 1e3, block.state[4] * 1e-2, block.state[5] * 1e-2,
          block.state[6] * 1e-2;
      r0_cel.block<3, 1>(0, 0) = t2c * r0_geo.block<3, 1>(0, 0);
      r0_cel.block<3, 1>(3, 0) = t2c * r0_geo.block<3, 1>(3, 0);
      //printf("Geocentric:"); for (int i=0; i<6; i++) printf(" %+12.4f ", r0_geo(i));
      //printf("\nInertial  :"); for (int i=0; i<6; i++) printf(" %+12.4f ", r0_cel(i));

      yPhi.block<6, 1>(0, 0) = r0_cel;

      // set the Phi part (aka, the state transition matrix) to I(6x6)
      int offset = 0;
      for (int col = 1; col < 7; col++) {
        int index = 6 * col + offset++;
        yPhi(index) = 1e0;
      }

      // t = current_time - t0 [seconds]
      double t = block.t.delta_sec(t0).to_fractional_seconds();
      // tout = t + step [seconds]
      double tout = t + step.to_fractional_seconds();

      if (Integrate) {
        // integrate
        Integrator.flag() = 1;
        Integrator.de(t, tout, yPhi, sol);

        // output epoch as datetime
        epoch = block.t;
        epoch.add_seconds(dso::milliseconds(static_cast<long>(tout * 1e3)));
        dso::strftime_ymd_hmfs(epoch, buf);

        printf("%s %+15.4f %+15.4f %+15.4f %+15.7f %+15.7f %+15.7f %15.5f\n",
               buf, sol(0), sol(1), sol(2), sol(3), sol(4), sol(5),
               epoch.as_mjd());
      } else {

        epoch = block.t;
        dso::strftime_ymd_hmfs(epoch, buf);
        printf("%s %+15.4f %+15.4f %+15.4f %+15.7f %+15.7f %+15.7f %15.5f\n",
               buf, r0_cel(0), r0_cel(1), r0_cel(2), r0_cel(3), r0_cel(4),
               r0_cel(5), epoch.as_mjd());
      }
    }
    ++rec_count;

    if (epoch.delta_sec(t0).to_fractional_seconds() > 60 * 60 * 12e0)
      break;

  } while (!error);

  printf("Num of records read: %6lu\n", rec_count);

  return (error == -1);
}
