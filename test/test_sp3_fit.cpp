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
#include "planetpos.hpp"
#include "sp3/sp3.hpp"
#include <cstdio>
#include <datetime/dtfund.hpp>

using dso::sp3::SatelliteId;

const int Integrate = true;
const int include_third_body = false;
// const int useVariationalEquations = false;

// approximate number of data points in a bulletin B file (disregarding
// preliminery values)
constexpr const int BBSZ = 40;

static unsigned call_nr = 0;

// usually using these datetimes ...
using Datetime = dso::datetime<dso::nanoseconds>;

// to transfer parameters for Variational Equations
struct AuxParams {
  double mjd_tai;
  double xp, yp, dut1;
  dso::HarmonicCoeffs *hc;
  dso::Mat2D<dso::MatrixStorageType::Trapezoid> *V, *W;
  int degree, order;
};

struct EopInfo {
  // actual size of arrays
  int sz;
  // arrays of EOP values extracted from C04
  double mjd[BBSZ], xpa[BBSZ], ypa[BBSZ], ut1a[BBSZ];
};

dso::Mat3x3 RzMat(double angle) noexcept {
  const double s = std::sin(angle);
  const double c = std::cos(angle);
  return dso::Mat3x3({c,s,0e0, -s,c,0e0, 0e0,0e0,1e0});
}

int getMeEops(const Datetime &t, char *bulletinb_fn, EopInfo *eopLUT,
              int download = 1) noexcept {
  if (download) {
    // get the corresponding bulletin b file from IERS (download)
    int error = dso::download_iers_bulletinb_for(t.mjd().as_underlying_type(),
                                                 bulletinb_fn);
    if (error) {
      fprintf(stderr, "Failed to download Bulletin B file for date: %ld\n",
              t.mjd().as_underlying_type());
      return 1;
    }
  }

  // parse the Bulletin B file (diregard preliminary values, only use final)
  dso::IersBulletinB_Section1Block bbblocks[35];
  dso::IersBulletinB bulb(bulletinb_fn);
  int bbblocks_size = bulb.parse_section1(bbblocks, false);
  if (bbblocks_size <= 0) {
    fprintf(stderr, "Failed parsing Bulletin B file %s\n", bulletinb_fn);
    return 2;
  }

  for (int i = 0; i < bbblocks_size; i++)
    eopLUT->mjd[i] = static_cast<double>(bbblocks[i].mjd);
  for (int i = 0; i < bbblocks_size; i++)
    eopLUT->xpa[i] = bbblocks[i].x;
  for (int i = 0; i < bbblocks_size; i++)
    eopLUT->ypa[i] = bbblocks[i].y;
  for (int i = 0; i < bbblocks_size; i++)
    eopLUT->ut1a[i] = bbblocks[i].dut1;

  eopLUT->sz = bbblocks_size;

  return 0;
}

// xp, yp in milliarcsecond [mas] dut1 in millisecond [ms]
int getEop(const Datetime &t, const EopInfo *eops, double &xp, double &yp,
           double &dut1) noexcept {
  // apply corrections; first use INTERP to interpolate x,y,ut1 values for
  // the given date (in [mas], [msec])
  if (iers2010::interp_pole(eops->mjd, eops->xpa, eops->ypa, eops->ut1a,
                            eops->sz, t.as_mjd(), xp, yp, dut1)) {
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

// handle gravity field
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

Eigen::Matrix<double, 3, 3>
ter2cel(double mjd_tai, double xp, double yp, double dut1,
        Eigen::Matrix<double, 3, 3> *dt2c) noexcept {

  // keep small part, do computations with this
  double mjd_days;
  const double taif = std::modf(mjd_tai, &mjd_days);
  const int leap_sec =
      dso::dat(dso::modified_julian_day(static_cast<int>(mjd_tai)));
  const double utcf = taif - static_cast<double>(leap_sec) / 86400e0;
  const double ut1f = utcf + (dut1 / 86400e0) * 1e-3;
  const double ttf  = taif + (32.184e0/86400e0);

  // const double mjd_utc = mjd_days + utcf;
  const double mjd_ut1 = mjd_days + ut1f;
  const double mjd_tt  = mjd_days + ttf;

  /*
  auto mat =
      iers2010::sofa::c2t06a(dso::mjd0_jd, mjd_tt, dso::mjd0_jd, mjd_ut1,
                             iers2010::DMAS2R * xp, iers2010::DMAS2R * yp);
  */
  // Form the celestial-to-intermediate matrix for this TT.
  auto rc2i = iers2010::sofa::c2i06a(dso::mjd0_jd, mjd_tt);
  // Predict the Earth rotation angle for this UT1.
  const double era = iers2010::sofa::era00(dso::mjd0_jd, mjd_ut1);
  // Estimate s'.
  const double sp = iers2010::sofa::sp00(dso::mjd0_jd, mjd_tt);
  // Form the polar motion matrix.
  auto rpom =
      iers2010::sofa::pom00(xp * iers2010::DMAS2R, yp * iers2010::DMAS2R, sp);
  // Combine to form the celestial-to-terrestrial matrix.
  auto mat = iers2010::sofa::c2tcio(rc2i, era, rpom);
  Eigen::Matrix<double, 3, 3> t2c(mat.data);
  
  /* ERA derivative */
  if (dt2c) {
    const dso::Mat3x3 S({0e0, iers2010::OmegaEarth, 0e0, -iers2010::OmegaEarth,
                         0e0, 0e0, 0e0, 0e0, 0e0});
    mat = rpom * (S*RzMat(era)) * rc2i;
    *dt2c = Eigen::Matrix<double, 3, 3>(mat.data);
  }

  return t2c;
}

void SunMoon(double mjd_tt, Eigen::Matrix<double, 3, 1> &sun,
             Eigen::Matrix<double, 3, 1> &moon) {
  const double jd = mjd_tt + dso::mjd0_jd; // date as jd (TT)
  double vec[3];
  dso::cspice::j2planet_pos_from(dso::cspice::jd2et(jd), 10, 399, vec);
  sun(0, 0) = vec[0] * 1e3;
  sun(1, 0) = vec[1] * 1e3;
  sun(2, 0) = vec[2] * 1e3;

  dso::cspice::j2planet_pos_from(dso::cspice::jd2et(jd), 301, 399, vec);
  moon(0, 0) = vec[0] * 1e3;
  moon(1, 0) = vec[1] * 1e3;
  moon(2, 0) = vec[2] * 1e3;

  return;
}

void VariationalEquations(
    double tsec,
    // state and state transition matrix
    const Eigen::VectorXd &yPhi,
    // state derivative and state transition matrix derivative
    Eigen::Ref<Eigen::VectorXd> yPhiP,
    //
    void *pAux) noexcept {

  // void pointer to AuxParameters
  AuxParams *params = static_cast<AuxParams *>(pAux);

  // current mjd
  const double cmjd = params->mjd_tai + tsec / dso::sec_per_day;
  // printf("\t> computing Variational Equations for t=%.6f, aka Mjd=%.6f call # %u\n", tsec, cmjd, call_nr);
  
  // terretrial to celestial for epoch
  Eigen::Matrix<double, 3, 3> t2c(
      ter2cel(cmjd, params->xp, params->yp, params->dut1, nullptr));
  // Eigen::Matrix<double, 3, 3> c2t = t2c.transpose();

  // split position and velocity vectors (inertial)
  Eigen::Matrix<double, 3, 1> r = yPhi.block<3, 1>(0, 0);
  Eigen::Matrix<double, 3, 1> v = yPhi.block<3, 1>(3, 0);

  // compute gravity-induced acceleration and gradient
  Eigen::Matrix<double, 3, 1> r_geo = t2c.transpose() * r;
  Eigen::Matrix<double, 3, 3> gpartials;
  Eigen::Matrix<double, 3, 1> gacc = dso::grav_potential_accel(
      r_geo, params->degree, params->order, *(params->V), *(params->W),
      *(params->hc), gpartials);

  // fucking crap! gravity acceleration are in earth-fixed frame; need to
  // have inertial acceleration!
  gacc = t2c * gacc;
  gpartials = t2c.transpose() * gpartials * t2c;

  // compute sun/moon induced acceleration
  Eigen::Matrix<double, 3, 3> tb_partials = Eigen::Matrix<double, 3, 3>::Zero();
  Eigen::Matrix<double, 3, 1> sacc = Eigen::Matrix<double, 3, 1>::Zero();
  Eigen::Matrix<double, 3, 1> macc = Eigen::Matrix<double, 3, 1>::Zero();
  if (include_third_body) {
    Eigen::Matrix<double, 3, 3> partials;
    Eigen::Matrix<double, 3, 1> rsun, rmoon;
    SunMoon(cmjd, rsun, rmoon);
    sacc =
        dso::point_mass_accel(params->hc->GM(), r, rsun, partials);
    tb_partials = partials;
    macc =
        dso::point_mass_accel(params->hc->GM(), r, rmoon, partials);
    tb_partials += partials;
  }

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
  dfdy.block<3, 3>(3, 0) = gpartials + tb_partials;
  dfdy.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Zero();

  // Derivative of combined state vector and state transition matrix
  // dPhi = dfdy * Phi;
  Eigen::Matrix<double, 6, 7> yPhip;
  yPhip.block<6, 6>(0, 1) = dfdy * Phi;

  // state derivative (aka [v,a]), in one (first) column
  yPhip.block<3, 1>(0, 0) = v;
  yPhip.block<3, 1>(3, 0) = gacc + sacc + macc;

  // matrix to vector (column-wise)
  yPhiP = Eigen::VectorXd(
      Eigen::Map<Eigen::VectorXd>(yPhip.data(), yPhip.cols() * yPhip.rows()));

  ++call_nr;
  return;
}

int main(int argc, char *argv[]) {

  // handle gravity field
  if (argc != 7) {
    fprintf(stderr,
            "Usage: %s <SP3 FILE> <GRAVITY MODEL FILE> [DEGREE] [ORDER] [SPK] "
            "[LSK]\n",
            argv[0]);
    return 1;
  }

  // parse degree and order of gravity field
  int degree = std::atoi(argv[3]); 
  int order = std::atoi(argv[4]);

  // to compute the planet's position via cspice, we need to load:
  // 1. the planetary ephemeris (SPK) kernel
  // 2. the leap-second (aka LSK) kernel
  dso::cspice::load_if_unloaded_spk(argv[5]);
  dso::cspice::load_if_unloaded_lsk(argv[6]);

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

  SatelliteId sv;
  dso::Sp3DataBlock block;

  if (sp3.num_sats() == 1) {
    sv.set_id(sp3.sattellite_vector()[0].id);
  } else {
    fprintf(stderr, "More than one Satellites in sp3 file. This test program "
                    "is only meant to work with one.\n");
    return 1;
  }

  // Handle EOPs
  printf("* downloading and parsing Bulletin B file ...\n");
  const auto t0 = sp3.start_epoch();
  char bulletinb_fn[256];
  EopInfo EopLookUpTables;
  if (getMeEops(t0, bulletinb_fn, &EopLookUpTables)) {
    fprintf(stderr, "Error. Failed to fetch EOPs data\n");
    return 3;
  }
  double xp, yp, dut1;
  if (getEop(t0, &EopLookUpTables, xp, yp, dut1)) {
    fprintf(stderr, "Error. Failed to interpolate EOPs\n");
    return 3;
  }
  printf("\t-> got eop data for date!\n");
  printf("\txp=%.15f yp=%.15f dut1=%.15f [mas], [msec]\n", xp, yp, dut1);

  // Create Auxiliary parameters
  AuxParams params{0e0, 0e0, 0e0, 0e0, &hc, &V, &W, degree, order};
  // SetUp an integrator
  printf("* setting up integrator ...\n");
  const double relerr = 1.0e-9  * 1e3;  // Relative and absolute
  const double abserr = 1.0e-16 * 1e2;  // accuracy requirement
  dso::SGOde Integrator(VariationalEquations, 6 * 7, relerr, abserr, &params);
  Eigen::Matrix<double, 6 * 7, 1> yPhi;
  Eigen::VectorXd sol(6 * 7);
  Eigen::VectorXd r0_geo(6), r0_cel(6);

  // Time
  // (remember) const auto t0 = sp3.start_epoch();
  auto epoch = t0; // current
  const auto step = sp3.interval();

  // let's try reading the records; note that -1 denotes EOF
  int error = 0;
  std::size_t rec_count = 0;
  char buf[64];
  printf("* iterating ...\n");
  do {

    // read next record ...
    if ((error = sp3.get_next_data_block(sv, block)) > 0) {
      printf("Something went wrong ....status = %3d\n", error);
      return 1;
    }

    // check the heath status
    bool position_ok = !block.flag.is_set(dso::Sp3Event::bad_abscent_position);

    if (position_ok) {
      
      // accumulate state (m, m/sec)
      r0_geo << block.state[0] * 1e3, block.state[1] * 1e3,
          block.state[2] * 1e3, block.state[4] * 1e-2, block.state[5] * 1e-2,
          block.state[6] * 1e-2;

      // transform geocentric state to inertial
      Eigen::Matrix<double, 3, 3> dt2c;
      Eigen::Matrix<double, 3, 3> t2c(
          ter2cel(block.t.as_mjd(), xp, yp, dut1, &dt2c));
      r0_cel.block<3, 1>(0, 0) = t2c * r0_geo.block<3, 1>(0, 0);
      r0_cel.block<3, 1>(3, 0) =
          t2c * r0_geo.block<3, 1>(3, 0) + dt2c * r0_geo.block<3, 1>(0, 0);

      // Vector containing state + variational equations
      yPhi = Eigen::Matrix<double, 6 * 7, 1>::Zero();
      yPhi.block<6, 1>(0, 0) = r0_cel;

      // set the Phi part (aka, the state transition matrix) to I(6x6)
      int offset = 0;
      for (int col = 1; col < 7; col++) {
        int index = 6 * col + offset++;
        yPhi(index) = 1e0;
      }

      params.mjd_tai = block.t.as_mjd();
      double tout = step.to_fractional_seconds();
      double t = 0e0;

      if (Integrate) {

        dso::strftime_ymd_hmfs(block.t, buf);
        
        // integrate (in inertial RF)
        Integrator.flag() = 1;
        Integrator.de(t, tout, yPhi, sol);
        if (std::abs(tout - t) > 5) {
          fprintf(stderr,
                  "Warning! Interpolation ended more than 5 secs away target: "
                  "%.6f reached: %.6f !\n",
                  tout, t);
        }

        // output epoch as datetime
        epoch = block.t;
        // THIS IS WRONG! the integrator has reached point t and NOT 
        // neccessarily tout
        // epoch.add_seconds(dso::milliseconds(static_cast<long>(tout * 1e3)));
        epoch.add_seconds( dso::nanoseconds(static_cast<long>(t * 1e9)) );
        dso::strftime_ymd_hmfs(epoch, buf);
        
        // inertial to terrestrial for SP3 comparisson
        t2c = ter2cel(epoch.as_mjd(), xp, yp, dut1, &dt2c);
        r0_geo.block<3, 1>(0, 0) = t2c.transpose() * sol.block<3, 1>(0, 0);
        r0_geo.block<3, 1>(3, 0) = t2c.transpose() * sol.block<3, 1>(3, 0) +
                                   dt2c.transpose() * sol.block<3, 1>(0, 0);

        // print results in terestrial coordinates
        printf("%s %+15.4f %+15.4f %+15.4f %+15.7f %+15.7f %+15.7f %18.9f\n",
               buf, r0_geo(0), r0_geo(1), r0_geo(2), r0_geo(3), r0_geo(4),
               r0_geo(5), epoch.as_mjd());
      } else {

        epoch = block.t;
        dso::strftime_ymd_hmfs(epoch, buf);

        printf("%s %+15.4f %+15.4f %+15.4f %+15.7f %+15.7f %+15.7f %18.9f\n",
               buf, r0_geo(0), r0_geo(1), r0_geo(2), r0_geo(3), r0_geo(4),
               r0_geo(5), epoch.as_mjd());
      }
    }
    ++rec_count;
    call_nr = 0;

    if (epoch.delta_sec(t0).to_fractional_seconds() > 12*3600e0)
      break;

  } while (!error && rec_count < 1000);

  printf("Num of records read: %6lu\n", rec_count);

  return (error == -1);
}
