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
#include <cassert>
#include <cstdio>
#include <datetime/dtfund.hpp>

using dso::sp3::SatelliteId;

const int Integrate = true;
const int include_third_body = true;
const int include_srp = false;

// approximate number of data points in a bulletin B file (disregarding
// preliminery values)
constexpr const int BBSZ = 40;

// count calls to VE
static unsigned call_nr = 0;

// Radiation pressure scale coefficient for H2C
// see https://ids-doris.org/documents/BC/satellites/DORISSatelliteModels.pdf
constexpr const double Cr_H2C = 0.88e0;
constexpr const double Mass_H2C = 1677e0; // [kg] total mass of H2C
constexpr const double Area_H2C = 39.71e0; // [m^2]

// Standard gravitational parameters for Sun and Moon in [km^3 / sec^2]
double GMSun, GMMoon;

// usually using these datetimes ...
using Datetime = dso::datetime<dso::nanoseconds>;

struct EopInfo {
  // actual size of arrays
  int sz;
  // arrays of EOP values extracted from C04
  double mjd[BBSZ], xpa[BBSZ], ypa[BBSZ], ut1a[BBSZ];
};

// to transfer parameters for Variational Equations
struct AuxParams {
  double mjd_tai;
  EopInfo *eopLookUpTables;
  dso::HarmonicCoeffs *hc;
  dso::Mat2D<dso::MatrixStorageType::Trapezoid> *V, *W;
  int degree, order;
};

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

// return xp, yp in milliarcsecond [mas]
// and dut1 in millisecond [ms]
int getEop(double mjd_utc, const EopInfo *eops, double &xp, double &yp,
           double &dut1) noexcept {
  // apply corrections; first use INTERP to interpolate x,y,ut1 values for
  // the given date (in [mas], [msec])
  if (iers2010::interp_pole(eops->mjd, eops->xpa, eops->ypa, eops->ut1a,
                            eops->sz, mjd_utc, xp, yp, dut1)) {
    fprintf(stderr, "ERROR. Failed call to interp_pole\n");
    return 1;
  }

  // account for variations in polar motion (Dx,Dy) ocean-tides; results in
  // [μas] and [μsec]
  double dxp, dyp, dut2;
  if (iers2010::ortho_eop(mjd_utc, dxp, dyp, dut2)) {
    fprintf(stderr, "ERROR. Failed call to ortho_eop!\n");
    return 4;
  }

  // transform to [mas] [msec]
  xp += dxp * 1e-3;
  yp += dyp * 1e-3;
  dut1 += dut2 * 1e-3;

  // account for libration effects; results in [μas]
  if (iers2010::pmsdnut2(mjd_utc, dxp, dyp)) {
    fprintf(stderr, "ERROR Failed call to pmsdnut2!\n");
    return 4;
  }

  // transform to [mas]
  xp += dxp * 1e-3;
  yp += dyp * 1e-3;

  return 0;
}

int getEop(const Datetime &t, const EopInfo *eops, double &xp, double &yp,
           double &dut1) noexcept {
  return getEop(t.as_mjd(), eops, xp, yp, dut1);
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
ter2cel(double mjd_tai, const EopInfo *eopLUT,
        Eigen::Matrix<double, 3, 3> *dt2c) noexcept {

  // keep small part, do computations with this
  double mjd_days;
  const double taif = std::modf(mjd_tai, &mjd_days);
  const int leap_sec =
      dso::dat(dso::modified_julian_day(static_cast<int>(mjd_tai)));
  const double utcf = taif - static_cast<double>(leap_sec) / 86400e0;

  // we now have date in UTC, get EOP values
  double xp, yp, dut1;
  if (getEop(mjd_days + utcf, eopLUT, xp, yp, dut1)) {
    fprintf(stderr, "ERROR. Failed getting EOP values");
  }

  const double ut1f = utcf + (dut1 / 86400e0) * 1e-3;
  const double ttf = taif + (32184e-3 / 86400e0);

  // const double mjd_utc = mjd_days + utcf;
  const double mjd_ut1 = mjd_days + ut1f;
  const double mjd_tt = mjd_days + ttf;

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
  // auto mat = iers2010::sofa::c2tcio(rc2i, era, rpom);
  auto mat = rpom * dso::Mat3x3::RotZ(era) * rc2i;
  
  // note that the following will result in an Eigen matrix that is the
  // transpose of mat (Eigen uses Column-Major and Mat3x3 Row-Major)
  Eigen::Matrix<double, 3, 3> t2c(mat.data);

  // ERA derivative
  if (dt2c) {
    const dso::Mat3x3 S({0e0, iers2010::OmegaEarth, 0e0, -iers2010::OmegaEarth,
                         0e0, 0e0, 0e0, 0e0, 0e0});
    mat = rpom * (S * dso::Mat3x3::RotZ(era)) * rc2i;
    *dt2c = Eigen::Matrix<double, 3, 3>(mat.data);
  }

  return t2c;
}

void SunMoon(double mjd_tai, const Eigen::Matrix<double, 3, 1> &rsat,
             Eigen::Matrix<double, 3, 1> &sun_acc,
             Eigen::Matrix<double, 3, 1> &moon_acc,
             Eigen::Matrix<double, 3, 1> &sun_pos,
             Eigen::Matrix<double, 3, 3> &mon_partials) noexcept {

  // TAI to TT (MJD)
  double mjd_days;
  const double taif = std::modf(mjd_tai, &mjd_days);
  const double ttf = taif + (32184e-3 / 86400e0);
  const double mjd_tt = mjd_days + ttf;

  const double jd = mjd_tt + dso::mjd0_jd; // date as JD (TT)
  double rsun[3], rmon[3];

  // position vector of sun/moon, in J2000, [km]
  dso::cspice::j2planet_pos_from(dso::cspice::jd2et(jd), 10, 399, rsun);
  dso::cspice::j2planet_pos_from(dso::cspice::jd2et(jd), 301, 399, rmon);

  Eigen::Matrix<double, 3, 1> rSun(rsun); // [km]
  Eigen::Matrix<double, 3, 1> rMon(rmon); // [km]

  // Sun-induced acceleration [km/sec^2]
  sun_acc = dso::point_mass_accel(GMSun, rsat * 1e-3, rSun);
  sun_acc = sun_acc * 1e-3; // [m/sec^2]

  // Moon-induced acceleration [m/sec^2]
  moon_acc = dso::point_mass_accel(GMMoon * 1e9, rsat, rMon * 1e3, mon_partials);

  // Sun position in [m]
  sun_pos = rSun * 1e3;

  return;
}

void VariationalEquations(
    double tsec, // TAI
    // state and state transition matrix (inertial RF)
    const Eigen::VectorXd &yPhi,
    // state derivative and state transition matrix derivative (inertial RF)
    Eigen::Ref<Eigen::VectorXd> yPhiP,
    // auxiliary parametrs
    void *pAux) noexcept {

  // void pointer to AuxParameters
  AuxParams *params = static_cast<AuxParams *>(pAux);

  // current mjd, TAI
  const double cmjd = params->mjd_tai + tsec / dso::sec_per_day;

  // terretrial to celestial for epoch
  Eigen::Matrix<double, 3, 3> t2c(
      ter2cel(cmjd, params->eopLookUpTables, nullptr));

  // split position and velocity vectors (inertial)
  Eigen::Matrix<double, 3, 1> r = yPhi.block<3, 1>(0, 0);
  Eigen::Matrix<double, 3, 1> v = yPhi.block<3, 1>(3, 0);

  // compute gravity-induced acceleration
  Eigen::Matrix<double, 3, 3> gpartials;
  Eigen::Matrix<double, 3, 1> r_geo = t2c.transpose() * r;
  Eigen::Matrix<double, 3, 1> gacc = dso::grav_potential_accel(
      r_geo, params->degree, params->order, *(params->V), *(params->W),
      *(params->hc), gpartials);

  // fucking crap! gravity acceleration in earth-fixed frame; need to
  // have inertial acceleration!
  gacc = t2c * gacc;
  gpartials = t2c.transpose() * gpartials * t2c;

  // third body perturbations, Sun and Moon [m/sec^2] in celestial RF
  Eigen::Matrix<double, 3, 1> rsun; // position of sun, [m] in celestial RF
  Eigen::Matrix<double, 3, 1> sun_acc;
  Eigen::Matrix<double, 3, 1> mon_acc;
  Eigen::Matrix<double, 3, 3> tb_partials;
  SunMoon(cmjd, r, sun_acc, mon_acc, rsun,
          tb_partials);
  if (!include_third_body) {
    sun_acc = Eigen::Matrix<double, 3, 1>::Zero();
    mon_acc = Eigen::Matrix<double, 3, 1>::Zero();
    tb_partials = Eigen::Matrix<double, 3, 3>::Zero();
  }

  // SRP
  Eigen::Matrix<double, 3, 1> srp = Eigen::Matrix<double, 3, 1>::Zero();
  if (include_srp) {
    dso::Vector3 rV({r(0), r(1), r(2)});
    dso::Vector3 sV({rsun(0), rsun(1), rsun(2)});
    if (utest::montebruck_shadow(rV, sV)) {
      srp = dso::solar_radiation_acceleration(r, rsun, Area_H2C, Mass_H2C,
                                              Cr_H2C);
    }
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
  yPhip.block<3, 1>(3, 0) = gacc + sun_acc + mon_acc + srp;

  // matrix to vector (column-wise)
  yPhiP = Eigen::VectorXd(
      Eigen::Map<Eigen::VectorXd>(yPhip.data(), yPhip.cols() * yPhip.rows()));

  //printf("# Acceleration components:\n");
  //printf("\t# Grav: %+.9f Sun: %+.9f Moon: %+.9f SRP: %+.9f [m/sec^2]\n",
  //       gacc(0), sun_acc(0), moon_acc(0), srp(0));
  //printf("\t# Grav: %+.9f Sun: %+.9f Moon: %+.9f SRP: %+.9f [m/sec^2]\n",
  //       gacc(1), sun_acc(1), moon_acc(1), srp(1));
  //printf("\t# Grav: %+.9f Sun: %+.9f Moon: %+.9f SRP: %+.9f [m/sec^2]\n",
  //       gacc(2), sun_acc(2), moon_acc(2), srp(2));

  ++call_nr;
  return;
}

int main(int argc, char *argv[]) {

  // handle gravity field
  if (argc < 8 || argc > 9) {
    fprintf(stderr,
            "Usage: %s <SP3 FILE> <GRAVITY MODEL FILE> [DEGREE] [ORDER] [SPK] [LSK] [PCK] [INTEGRATION INTERVAL IN SEC - optional-]"
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
  // 3. the planetary constants kernel (PCK)
  dso::cspice::load_if_unloaded_spk(argv[5]);
  dso::cspice::load_if_unloaded_lsk(argv[6]);
  // well, actually we do not need to load the kernel; just get the values we 
  // want!
  // dso::cspice::load_if_unloaded_pck(argv[7]);
  dso::getSunMoonGM(argv[7], GMSun, GMMoon); // [km^3 / sec^2]

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

  // a satellite instance, will hold the sp3 satellite
  SatelliteId sv;
  dso::Sp3DataBlock block;

  if (sp3.num_sats() == 1) {
    sv.set_id(sp3.sattellite_vector()[0].id);
  } else {
    fprintf(stderr, "More than one Satellites in sp3 file. This test program "
                    "is only meant to work with one.\n");
    return 1;
  }

  // Handle EOPs:
  // Download bulletin B for the date, parse file and store Look Up Tables
  // in an EopInfo instance (for xp, yp, dut1)
  printf("* downloading and parsing Bulletin B file ...\n");
  const auto t0 = sp3.start_epoch();
  char bulletinb_fn[256];
  EopInfo EopLookUpTables;
  if (getMeEops(t0, bulletinb_fn, &EopLookUpTables)) {
    fprintf(stderr, "Error. Failed to fetch EOPs data\n");
    return 3;
  }

  // Create Auxiliary parameters. This will be passed-in the variational
  // equations to compute derivs.
  AuxParams params{
      sp3.start_epoch().as_mjd(), &EopLookUpTables, &hc, &V, &W, degree, order};

  // SetUp an integrator
  printf("* setting up integrator ...\n");
  const double relerr = 1e-12;  // Relative and absolute
  const double abserr = 1e-12; // accuracy requirement
  dso::SGOde Integrator(VariationalEquations, 6 + 6*6, relerr, abserr, &params);

  // matrices ...
  Eigen::Matrix<double, 6+6*6, 1> yPhi;     // state + var. equations
  Eigen::VectorXd sol(6+6*6);               // integrator solution
  Eigen::Matrix<double, 6*6,1> I6x6_vec;
  Eigen::VectorXd r0_geo(6), r0_cel(6);

  // Time
  // (remember) const auto t0 = sp3.start_epoch();
  auto epoch = t0;                  // current
  auto step = sp3.interval(); // integration interval
  if (argc == 8) {
      step = dso::nanoseconds(std::atoi(argv[7]) * 1'000'000'000L);
  }

  // let's try reading the records; note that -1 denotes EOF
  int error = 0;
  std::size_t rec_count = 0;
  char buf[64];
  printf("* iterating ...\n");
  do {

    // read next record ...
    if ((error = sp3.get_next_data_block(sv, block)) > 0) {
      fprintf(stderr, "Something went wrong ....status = %3d\n", error);
      return 1;
    }

    // check the heath status
    bool position_ok = !block.flag.is_set(dso::Sp3Event::bad_abscent_position);

    if (position_ok) {

      // accumulate state (m, m/sec), earth-fixed
      r0_geo << block.state[0] * 1e3, block.state[1] * 1e3,
          block.state[2] * 1e3, block.state[4] * 1e-1, block.state[5] * 1e-1,
          block.state[6] * 1e-1;

      // Terrestrial to Celestial transformation matrix and derivative
      Eigen::Matrix<double, 3, 3> dt2c;
      Eigen::Matrix<double, 3, 3> t2c(
          ter2cel(block.t.as_mjd(), params.eopLookUpTables, &dt2c));

      // transform geocentric state to inertial
      r0_cel.block<3, 1>(0, 0) = t2c * r0_geo.block<3, 1>(0, 0);
      r0_cel.block<3, 1>(3, 0) =
          t2c * r0_geo.block<3, 1>(3, 0) + dt2c * r0_geo.block<3, 1>(0, 0);

      if (Integrate) {

        // Vector containing state + variational equations size: 6 + 6x6
        // Ref. Frame: inertial
        yPhi = Eigen::Matrix<double, 6 + 6*6, 1>::Zero();
        yPhi.block<6,1>(0,0) = r0_cel;
        {
          int k = 6;
          for (int col=1; col<7; col++) for (int row=0; row<6; row++) yPhi(k++) = (col-1==row) ? 1e0 : 0e0;
        }

        // t0 for variational equations (TAI)
        params.mjd_tai = block.t.as_mjd();

        // target t for variational equations; seconds after t0
        double tout = step.to_fractional_seconds();

        // seconds after t0 for integrator
        double t = 0e0;

        // integrate (in inertial RF), from 0 to step
        // the "real" date (in mjd) inside the integrator is:
        // params.mjd_tai + t / 86400
        Integrator.flag() = 1;
        Integrator.de(t, tout, yPhi, sol);
        if (std::abs(tout - t) > 5) {
          fprintf(stderr,
                  "Warning! Interpolation ended more than 5 secs away target: "
                  "%.6f reached: %.6f\n",
                  tout, t);
        }
        // solution vector stored in sol
        // time is t seconds after initial t (which is 0), hence time
        // reached is: block.t(=params.mjd_tai) + t/86400

        // output epoch as datetime
        epoch = block.t;
        epoch.add_seconds(
            dso::nanoseconds(static_cast<unsigned long>(t * 1e9)));
        dso::strftime_ymd_hmfs(epoch, buf);

        // transform solution to terrestrial for SP3 comparisson
        t2c = ter2cel(epoch.as_mjd(), params.eopLookUpTables, &dt2c);
        r0_geo.block<3, 1>(0, 0) = t2c.transpose() * sol.block<3, 1>(0, 0);
        r0_geo.block<3, 1>(3, 0) = t2c.transpose() * sol.block<3, 1>(3, 0) +
                                   dt2c.transpose() * sol.block<3, 1>(0, 0);

        // print results in terestrial coordinates
        printf("%s %+15.4f %+15.4f %+15.4f %+15.7f %+15.7f %+15.7f %18.9f\n",
               buf, r0_geo(0), r0_geo(1), r0_geo(2), r0_geo(3), r0_geo(4),
               r0_geo(5), epoch.as_mjd());

        // print results in celestial coordinates
        // printf("%s %+15.4f %+15.4f %+15.4f %+15.7f %+15.7f %+15.7f %18.9f\n",
        //       buf, sol(0), sol(1), sol(2), sol(3), sol(4),
        //       sol(5), epoch.as_mjd());

      } else {

        // Do not integrate, just report results for this epoch
        // Ref. Frame: Earth-fixed (as in sp3 file)
        // Time Scale: the one recorded in the sp3 file, probably TAI
        // Units     : [m] and [m/sec]

        epoch = block.t;
        dso::strftime_ymd_hmfs(epoch, buf);

        // print terrestrial satellite oordinates
        printf("%s %+15.4f %+15.4f %+15.4f %+15.7f %+15.7f %+15.7f %18.9f\n",
               buf, r0_geo(0), r0_geo(1), r0_geo(2), r0_geo(3), r0_geo(4),
               r0_geo(5), epoch.as_mjd());

        // print celestial satellite oordinates
        // printf("%s %+15.4f %+15.4f %+15.4f %+15.7f %+15.7f %+15.7f %18.9f\n",
        //       buf, r0_cel(0), r0_cel(1), r0_cel(2), r0_cel(3), r0_cel(4),
        //       r0_cel(5), epoch.as_mjd());
      }

      ++rec_count;
    }

    if (epoch.delta_sec(t0).to_fractional_seconds() > 12 * 3600e0)
      break;

  } while (!error && rec_count < 1000);

  printf("Num of records read: %6lu\n", rec_count);

  return (error == -1);
}
