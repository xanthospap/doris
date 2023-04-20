#include "astrodynamics.hpp"
#include "beacon_tbl.hpp"
#include "datetime/datetime_write.hpp"
#include "datetime/utcdates.hpp"
#include "doris_observation_equations.hpp"
#include "doris_rinex.hpp"
#include "doris_utils.hpp"
#include "filters/filters.hpp"
#include "geodesy/geoconst.hpp"
#include "geodesy/units.hpp"
#include "iers2010/cel2ter.hpp"
#include "iers2010/hardisp.hpp"
#include "iers2010/iers2010.hpp"
#include "iers2010/iersc.hpp"
#include "iers2010/tropo.hpp"
#include "integrators.hpp"
#include "orbit_integration.hpp"
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

constexpr const int INCLUDE_ATTITUDE = true;
constexpr const int INCLUDE_ATM_DRAG = true;
constexpr const int INCLUDE_SRP_DRAG = true;
constexpr const int m = INCLUDE_ATM_DRAG + INCLUDE_SRP_DRAG;
constexpr const int n = 6 + m;
/* only compute Doppler count if two observation are within this time interval
 */
constexpr const double RESTART_AFTER_SEC = 11e0;
constexpr const double MAX_HOURS = 6;
/* signal a new satellite pass over a beacon */
constexpr const double NEW_PASS_AFTER_MIN = 30e0;

int site_affected_by_ssa(const char *b4id) noexcept {
  int affected = 0;
  affected += !(std::strncmp("ARFB", b4id, 4));
  affected += !(std::strncmp("CADB", b4id, 4));
  return affected;
}

struct {
  const char *id;
  const char *fallback;
} FallBackVmfNames[] = {
  {"CRQC", "CRQB"},
  {"KEYC", "KEWC"},
  {"MALC", "MALB"},
  {"ROZC", "ROXC"},
  {"SJVC", "SJUC"},
  {"REVC", "REUB"},
  {"SVBC", "SVAC"}
};
constexpr const int FallBackVmfNamesSize =
    sizeof(FallBackVmfNames) / sizeof(FallBackVmfNames[0]);
const char *fallback_name(const char *fcid) noexcept {
  for (int i=0; i<FallBackVmfNamesSize; i++)
    if (!std::strncmp(fcid, FallBackVmfNames[i].id, 4))
      return FallBackVmfNames[i].fallback;
  return nullptr;
}

struct Kalman {
  dso::TwoPartDate t;
  Eigen::VectorXd x;
  Eigen::MatrixXd P;
  Eigen::VectorXd K;

  Kalman(int numParams, dso::TwoPartDate tai = dso::TwoPartDate{}) noexcept
      : t(tai), x(Eigen::VectorXd::Zero(numParams)),
        P(Eigen::MatrixXd::Identity(numParams, numParams)),
        K(Eigen::VectorXd::Zero(numParams)) {}

  void reset_frequency_bias(int index, double val0, double sigma0) noexcept {
    x(index) = val0;
    P.row(index).setZero();
    P.col(index).setZero();
    P(index,index) = sigma0*sigma0;
  }

  void time_update(const dso::TwoPartDate &ti,
                   const Eigen::MatrixXd &F) noexcept {
    // assert(F.rows() == F.cols() == P.rows());
    P = F * P * F.transpose();
    t = ti;
  }

  Eigen::VectorXd &estimates() noexcept { return x; }
  Eigen::VectorXd estimates() const noexcept { return x; }

  double parameter(int index) const noexcept { return x(index); }
  double &parameter(int index) noexcept { return x(index); }

  double prediction_residual(double z, double g, double sigma,
                             const Eigen::VectorXd &H, double &var) noexcept {
    const double R = sigma * sigma;
    var = R + H.dot(P * H);
    return (z - g) - H.transpose() * (K * (z - g));
  }
  
  auto observation_no_update(double z, double g, double sigma,
                          const Eigen::VectorXd &H) noexcept {
    double inv_w = sigma * sigma;
    // kalman gain
    K = P * H / (inv_w + H.dot(P * H));
    // state update (note that y is a scalar y <- z - g = Y(k) - G(X,t)
    return x + K * (z - g);
  }

  void observation_update(double z, double g, double sigma,
                          const Eigen::VectorXd &H) noexcept {
    double inv_w = sigma * sigma;
    // kalman gain
    K = P * H / (inv_w + H.dot(P * H));
    // state update (note that y is a scalar y <- z - g = Y(k) - G(X,t)
    x = x + K * (z - g);
    const int N = x.rows();
    // covariance update (Joseph variant)
    auto KWm1Kt = (K * sigma) * (K * sigma).transpose();
    auto ImKG = Eigen::MatrixXd::Identity(N, N) - K * H.transpose();
    P = ImKG * P * ImKG.transpose() + KWm1Kt;
  }
};

struct TropoDetails {
  // Dtropo_hydrostatic = zhd0 * mfh;
  // Dtropo_wet         = zwd0 * mfw;
  double zhd, mfh, zwd, mfw;
  double operator()() const noexcept { return zhd * mfh + zwd * mfw; }
};

class RcSv {
  dso::TwoPartDate tai{dso::TwoPartDate()}; /* time of reception */
  dso::TwoPartDate etai{dso::TwoPartDate()}; /* time of emission */
  dso::TwoPartDate tprop{dso::TwoPartDate()}; /* proper time (receiver) */
  char fcid[5] = {'\0'};
  double l2_cycles; /* L [2GHz] measurement */
  double dIon{0};   /* Iono correction in L2GHz cycles */
  TropoDetails dTro;
  Eigen::Matrix<double, 3, 1> svGcrf; /* Sv position at t=tai, GCRF [m] */
  Eigen::Matrix<double, 3, 1> bcItrf; /* Beacon coordinates, ITRF [m] */
  Eigen::Matrix<double, 3, 3> R; /* ITRF-to-GCRF matrix at t=etai */
  int _restart = {0};

public:
  // geometric distance:
  // d = | svGcrf - R(etai) * bcItrf |
  double s() const noexcept {
    return (svGcrf - R*bcItrf).norm();
  }

  void mark_restart(int onoff=1) noexcept { _restart = onoff; }

  dso::TwoPartDate time_of_reception() const noexcept {return tai;}
  dso::TwoPartDate proper_time() const noexcept {return tprop;}
  Eigen::Matrix<double, 3, 1> sv_gcrf() const noexcept {return svGcrf;}
  Eigen::Matrix<double, 3, 1> bc_itrf() const noexcept {return bcItrf;}
  TropoDetails tropo() const noexcept {return dTro;}
  double l2GHz_cycles() const noexcept {return l2_cycles;}
  double iono_cycles() const noexcept {return dIon;}
  const char *bc_name() const noexcept {return fcid;}
  int restart() const noexcept {return _restart;}

  void set_emission_tai(const dso::TwoPartDate &t,
                        const dso::EopLookUpTable &eops) noexcept {
    dso::Itrs2Gcrs Rot(t.tai2tt(), &eops);
    R = Rot.itrf2gcrf();
  }
  void set_proper_time(const dso::TwoPartDate &t) noexcept {
    tprop = t;
  }

  RcSv(const char *b4cid, const dso::TwoPartDate &t, double L2, double dion,
       const TropoDetails &trp, const Eigen::Matrix<double, 3, 1> &sv_gcrf,
       const Eigen::Matrix<double, 3, 1> &bc_itrf,
       const dso::EopLookUpTable &eops) noexcept
      : tai(t), etai(t), l2_cycles(L2), dIon(dion), dTro(trp), svGcrf(sv_gcrf),
        bcItrf(bc_itrf), _restart(0) {
    std::memcpy(fcid, b4cid, sizeof(char) * 4);
    dso::Itrs2Gcrs Rot(etai.tai2tt(), &eops);
    R = Rot.itrf2gcrf();
  }
};

dso::TwoPartDate correction_aberration(const RcSv &v, const dso::Itrs2Gcrs &R) noexcept {
  // all work done in GCRF
  assert(R.tt() == v.time_of_reception().tai2tt());
  dso::TwoPartDate tai_emission(v.time_of_reception());
  // time correction (sec)
  const double sec = v.s() / iers2010::C;
  // time of emission
  tai_emission._small -= (sec / 86400e0);
  tai_emission.normalize();
  return tai_emission;
}

/* compute Ue + Ve^2 / 2 */
double emitter_potential(const Eigen::Matrix<double, 3, 1> &ritrf,
                         double omegaEarth, double GMEarth) noexcept {
  const double p2 = ritrf(0) * ritrf(0) + ritrf(1) * ritrf(1);
  const double Ve22 = (omegaEarth * omegaEarth * p2) / 2e0;
  const double Ue = GMEarth / ritrf.norm();
  return Ue + Ve22;
}

int observation_equation(const RcSv &v1, const RcSv &v2, double feN, double frT,
                         double DfeFen, double emitterPotential, double &Vobs, double &Vtheo,
                         double &dDfeNfeN, Eigen::Matrix<double, 3, 1> &dObsdr,
                         Eigen::Matrix<double, 3, 1> &dObsdv) noexcept {
  constexpr const double c = iers2010::C;
  const double dt =
      v2.proper_time().diff<dso::DateTimeDifferenceType::FractionalSeconds>(
          v1.proper_time());
  const double dIon = (c / feN) * (v2.iono_cycles() - v1.iono_cycles()) / dt;
  const double dTro = (v2.tropo()() - v1.tropo()()) / dt;
  const double Ndop = v2.l2GHz_cycles() - v1.l2GHz_cycles();
  const double fac = 1e0 - emitterPotential/c/c;

  Vobs = (c / feN) * (feN - frT - Ndop / dt) + dIon;
  Vtheo =
      fac*(v2.s() - v1.s()) / dt - (c / feN) * (Ndop / dt + frT) * DfeFen + dTro;

  // partials
  dDfeNfeN = (c / feN) * (Ndop / dt + frT);
  // grad of Vtheo w.r.t dr (satellite in GCRF)
  { dObsdr = (v2.sv_gcrf() / v2.s() - v1.sv_gcrf() / v1.s()) / dt; }
  // grad of Vtheo w.r.t dv (satellite in GCRF)
  { dObsdv = Eigen::Matrix<double, 3, 1>::Zero(); }

  return 0;
}

int get_tropo_vmf(const char *fcid, const dso::datetime<dso::nanoseconds> &t,
                  const Eigen::Matrix<double, 3, 1> &bxyz, double zd,
                  dso::SiteVMF3Feed &feed, TropoDetails &Dtrop) noexcept {
  assert(zd >= 0e0 && zd <= dso::DPI / 2e0);

  // ellipsoidal coordinates of the station; store them in an array
  // note: ellipsoidal = longitude, latitude, h_ell
  Eigen::Matrix<double, 3, 1> bell = dso::car2ell<dso::ellipsoid::grs80>(bxyz);

  // interpolate for site
  dso::vmf3_details::SiteVMF3GRMeteoRecord meteo;
  if (feed.interpolate(fcid, t, meteo)) {
    // before quiting, try for a fallback site name
    const char *fallback_site_name = fallback_name(fcid);
    if (fallback_site_name) {
      printf("[DEBUG] Failed Tropo; trying for fallback name %.4s (%.4s) for VMF3 grid ...",
             fallback_site_name, fcid);
      if (feed.interpolate(fallback_site_name, t, meteo)) {
        printf(" Noop!\n");
        fprintf(stderr, "[WRNNG] Observation skipped for site %.4s (TROPO missing)\n", fcid);
        return 1;
      } else {
        printf(" yeap!\n");
      }
    } else {
      fprintf(stderr, "[WRNNG] Observation skipped for site %.4s (TROPO missing)\n", fcid);
      return 1;
    }
  }

  double mfh, mfw;
  if (dso::vmf3(meteo.ah, meteo.aw, t, bell(1), bell(0), dso::deg2rad(zd), mfh,
                mfw)) {
    fprintf(stderr, "ERROR Failed to compute vmf3!\n");
    return 1;
  }

  Dtrop.zhd = meteo.zhd;
  Dtrop.mfh = mfh;
  Dtrop.zwd = meteo.zwd;
  Dtrop.mfw = mfw;

  return 0;
}

std::vector<const char *>
beaconcrs2cchar_vec(const std::vector<dso::BeaconCoordinates> &bcv) noexcept {
  std::vector<const char *> site_names;
  site_names.reserve(bcv.size());
  for (auto it = bcv.begin(); it != bcv.end(); ++it) {
    site_names.push_back(it->id);
  }
  /* check if we have fallback site names and replace */
  for (auto it = site_names.begin(); it != site_names.end(); ++it) {
    for (int i=0; i<FallBackVmfNamesSize; i++) {
      if (!std::strncmp(*it, FallBackVmfNames[i].id, 4)) {
        printf("[NOTE] Replacing site name %.4s with %.4s in VMF3 feed\n", *it, FallBackVmfNames[i].fallback);
        *it = FallBackVmfNames[i].fallback;
      }
    }
  }
  return site_names;
}

auto findLastObs(const char *fcid, std::vector<RcSv> &vec) noexcept {
  return std::find_if(vec.begin(), vec.end(), [=](const RcSv &it) {
    return !(std::strncmp(fcid, it.bc_name(), 4));
  });
}
auto findItrfCrd(const char *fcid,
                 const std::vector<dso::BeaconCoordinates> &vec) noexcept {
  return std::find_if(vec.begin(), vec.end(),
                      [=](const dso::BeaconCoordinates &it) {
                        return !(std::strncmp(fcid, it.id, 4));
                      });
}
auto getItrfCrd(
    std::vector<dso::BeaconCoordinates>::const_iterator it) noexcept {
  return Eigen::Matrix<double, 3, 1>({it->x, it->y, it->z});
}

double iono_l2_correction(const dso::BeaconObservations &obs, int l1i,
                          int l2i) noexcept {
  /* (L[2GHz] - sqrt(γ) L[400MHz]) / (γ-1) in [cycles]*/
  return (obs.m_values[l1i].m_value -
          dso::GAMMA_FACTOR_SQRT * obs.m_values[l2i].m_value) /
         (dso::GAMMA_FACTOR - 1e0);
}

int check_obs_flags(const dso::BeaconObservations &obs, int l1i, int l2i,
                    int w1i, int w2i) noexcept {
  if (obs.m_values[l1i].m_flag1 == '1' || // 2 GHz central frequency measurement
      obs.m_values[l1i].m_flag2 == '1' || // discontinuity of 2 GHZ measurement
      obs.m_values[l2i].m_flag1 ==
          '1' || // 400 MHz central frequency measurement
      obs.m_values[l2i].m_flag2 ==
          '1' || // discontinuity of 400 MHZ measurement
      obs.m_values[w1i].m_flag1 == '1' || // station on restart mode
      obs.m_values[w2i].m_flag1 == '1'    // station on restart mode
  )
    return 1;
  return 0;
}

int get_rinex_indexes(const dso::DorisObsRinex &rnx, int &l1_idx, int &l2_idx,
                      int &f_idx, int &w1_idx, int &w2_idx) noexcept {
  l1_idx = rnx.get_observation_code_index(
      dso::ObservationCode{dso::ObservationType::phase, 1});

  // index of the 400MHz phase measurement (need for iono-free reduction)
  l2_idx = rnx.get_observation_code_index(
      dso::ObservationCode{dso::ObservationType::phase, 2});

  // index of the F measurement (relative frequency offset)
  f_idx = rnx.get_observation_code_index(
      dso::ObservationCode{dso::ObservationType::frequency_offset});

  // index of the W measurement (power level), frequency 2GHz
  w1_idx = rnx.get_observation_code_index(
      dso::ObservationCode{dso::ObservationType::power_level, 1});

  // index of the W measurement (power level), frequency 400MHz
  w2_idx = rnx.get_observation_code_index(
      dso::ObservationCode{dso::ObservationType::power_level, 2});

  if (f_idx < 0 || l1_idx < 0 || l2_idx < 0 || w1_idx < 0 || w2_idx < 0) {
    fprintf(stderr,
            "[ERROR] Failed to find requested Observation Types in RINEX\'s "
            "observation "
            "types vector! (traceback: %s)\n",
            __func__);
    return 1;
  }

  return 0;
}

// hold satellite state & time
class OrbitIntegrator {
  // Current datetime in TAI
  dso::TwoPartDate mjd_tai;
  // state vector at t=tai in GCRF, at DORIS receiver RP (Iono-Free)
  Eigen::Matrix<double, 6, 1> state;
  // state transition matrix at t=tai
  Eigen::Matrix<double, 6, 6> Phi;
  // an instance to transform between ITRF and GCRF coordinates
  dso::Itrs2Gcrs Rot;
  // SV details
  dso::SvFrame *sv_frame{nullptr};

public:
  OrbitIntegrator(const dso::TwoPartDate &tai, const dso::EopLookUpTable &eops)
      : mjd_tai(tai), Rot(tai.tai2tt(), &eops){};

  void set_attitude(const dso::IntegrationParameters &params) {
    sv_frame = params.svFrame;
  }

  dso::TwoPartDate &tai_time() noexcept { return mjd_tai; }
  dso::TwoPartDate tai_time() const noexcept { return mjd_tai; }

  Eigen::Matrix<double, 6, 6> stateTransitionMatrix() const noexcept {
    return Phi;
  }

  Eigen::Matrix<double, 3, 1> gcrf_position_cm() const noexcept {
    return state.block<3, 1>(0, 0);
  }
  Eigen::Matrix<double, 6, 1> gcrf_state_cm() const noexcept { return state; }
  Eigen::Matrix<double, 6, 1> &gcrf_state_cm() noexcept { return state; }
  Eigen::Matrix<double, 3, 1> itrf_position_cm() const noexcept {
    return Rot.gcrf2itrf(gcrf_position_cm());
  }
  Eigen::Matrix<double, 6, 1> itrf_state_cm() const noexcept {
    return Rot.gcrf2itrf(gcrf_state_cm());
  }
  
  Eigen::Matrix<double, 3, 1> gcrf_position_arp() noexcept {
    // return satellite's Antenna L2 RP position, in GCRF, assuming that 
    // state is referred to CoM (GCRF)
    if (sv_frame) {
      // GCRF: from satellite ARP to CoG
      Eigen::Matrix<double, 3, 1> dr;
      assert(! (sv_frame->arp2cog(tai_time(), dr)) );
      return gcrf_position_cm() - dr;
    }
    return gcrf_position_cm();
  }
  
  Eigen::Matrix<double, 3, 1> itrf_position_arp() noexcept {
    return Rot.gcrf2itrf(gcrf_position_arp());
  }

  /* assumes state is w.r.t SV's CoM */
  void set_state_from_itrf(const dso::TwoPartDate &tai,
                           const Eigen::Matrix<double, 6, 1> &itrf) noexcept {
    mjd_tai = tai;
    Rot.prepare(mjd_tai.tai2tt());
    gcrf_state_cm() = Rot.itrf2gcrf(itrf);
  }

  Eigen::Matrix<double, 6, 6> extractStateTransitionMatrix(
      const Eigen::Matrix<double, (6 + 6 * n), 1> &yP) noexcept {
    Eigen::Matrix<double, 6, 6> F;
    int of = 6;
    F.row(0) = yP.block<6, 1>(of, 0);
    F.row(1) = yP.block<6, 1>(of + 1 * n, 0);
    F.row(2) = yP.block<6, 1>(of + 2 * n, 0);
    F.row(3) = yP.block<6, 1>(of + 3 * n, 0);
    F.row(4) = yP.block<6, 1>(of + 4 * n, 0);
    F.row(5) = yP.block<6, 1>(of + 5 * n, 0);
    return F;
  }

  void
  setInitialconditions(Eigen::Matrix<double, (6 + 6 * n), 1> &yP0) noexcept {
    int of = 6;
    yP0.block<6 * n, 1>(6, 0) = Eigen::Matrix<double, (6 * n), 1>::Zero();
    for (int i = 0; i < 6; i++) {
      int start = of + i * n;
      yP0(start + i) = 1e0;
    }
  }

  int integrate(const dso::TwoPartDate &mjd_target,
                dso::SGOde &integrator) noexcept {

    // number of equations:
    // 6 for the state + 6*n for the variational equations, where n = 6 + m
    constexpr const int NumEqn = 6 + 6 * n;

    // Vector containing state + variational equations size: 6 + 6xn
    Eigen::Matrix<double, NumEqn, 1> yPhi0 =
        Eigen::Matrix<double, NumEqn, 1>::Zero();

    // transform state from ITRF to GCRF
    yPhi0.block<6, 1>(0, 0) = gcrf_state_cm();

    // Initial condition for state transition matrix Φ(t0,t0) = I
    setInitialconditions(yPhi0);

    // t0 for integration (TAI), aka current reference time
    integrator.params->reference_epoch() = mjd_tai.normalized();

    // target t for integration: seconds after t0
    double tout =
        mjd_target.diff<dso::DateTimeDifferenceType::FractionalSeconds>(
            integrator.params->reference_epoch());

    // integration solution
    Eigen::VectorXd sol(NumEqn);

    /* char buf[64];
    {
      using ST = dso::nanoseconds::underlying_type;
      const auto T = mjd_tai.normalized();
      dso::datetime<dso::nanoseconds> wt(
          dso::modified_julian_day(T._big),
          dso::nanoseconds(static_cast<ST>(
              (T._small * 86400e0) * dso::nanoseconds::sec_factor<double>())));
      dso::strftime_ymd_hmfs(wt, buf);
      fprintf(stderr, "[LOG] Integrating from %s to ", buf);
      wt.add_seconds(dso::nanoseconds(static_cast<ST>(
          (tout) * dso::nanoseconds::sec_factor<double>())));
      dso::strftime_ymd_hmfs(wt, buf);
      fprintf(stderr, "%s (%.9f sec away)\n", buf, tout);
    }*/

    // integrate (in inertial RF), from 0+mjd_tai to tout+mjd_tai [sec]
    double tsec = 0e0;
    integrator.flag() = dso::SGOde::IFLAG::RESTART;
    if (integrator.de(tsec, tout, yPhi0, sol) != dso::SGOde::IFLAG::SUCCESS) {
      fprintf(stderr, "[ERROR] Integrator error!\n");
      return 1;
    }

    // after integration, we reached the epoch: t0+tsec[sec] (TAI)
    const dso::TwoPartDate tout_mjd(
        dso::TwoPartDate(mjd_tai + dso::TwoPartDate(0e0, tsec / 86400e0))
            .normalized());

    // let's see were we are at
    if (std::abs(tout_mjd.diff<dso::DateTimeDifferenceType::FractionalSeconds>(
            mjd_target)) > 1e-9) {
      fprintf(stderr,
              "ERROR wanted integration to %.9f and got up to %.9f, that is "
              "%.9f sec apart!\n",
              mjd_target.mjd(), tout_mjd.mjd(),
              (mjd_target.mjd() - mjd_tai.mjd()) * 86400e0);
      return 1;
    }

    // everything seems ok, update state and time
    mjd_tai = tout_mjd;

    // assign Phi matrix 6x6
    Phi = extractStateTransitionMatrix(sol);
    state = sol.block<6, 1>(0, 0);

    // ! warning !
    // do not forget this step, rotation matrix for now
    Rot.prepare(mjd_tai.tai2tt());

    return 0;
  }
};

int main(int argc, char *argv[]) {
  // check input
  if (argc != 2) {
    fprintf(stderr, "USAGE: %s [YAML CONFIG]\n", argv[0]);
    return 1;
  }

  // resolve the yaml config file and get the root node
  const YAML::Node config = YAML::LoadFile(argv[1]);
  char buf[256];
  int error;

  // Elevation Cut-Off Angle in [degrees]
  // -------------------------------------------------------------------------
  double EleCutOff; // [deg]
  if (dso::get_yaml_value_depth2<double>(config, "filtering",
                                         "elevation-cut-off", EleCutOff)) {
    fprintf(stderr, "ERROR failed to find spk kernel\n");
    return 1;
    printf("[EDBUG] Using elevation cut-off angle = %.3f\n", EleCutOff);
  }

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

  // DORIS RINEX instance
  // -------------------------------------------------------------------------
  if (dso::get_yaml_value_depth2(config, "data", "doris-rinex", buf)) {
    fprintf(stderr, "ERROR. Failed parsing data/rinex file from YAML %s\n",
            argv[1]);
    return 1;
  }
  // Construct the Doris RINEX instance rnx
  dso::DorisObsRinex rnx(buf);
  // Fit a linear Model for the Relative Receiver Offset Values (use proper
  // time for this).
  dso::PolynomialModel<dso::datetime<dso::nanoseconds>> rfo_fit(1);
  {
    if (dso::fit_relative_frequency_offset(buf, rfo_fit, false, 4)) {
      fprintf(stderr, "ERROR! Failed to fit RFO values!\n");
      return 2;
    }
  }

  // Initial Orbit
  // -------------------------------------------------------------------------
  // Initial satellite state, get it from the Sp3 using the RINEX's time of
  // first obs
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
    const int ref_mjd = rnx.ref_datetime().as_mjd();
    const int start = ref_mjd - 8;
    const int end = ref_mjd + 10;
    // parse C04 EOPs
    if (parse_iers_C0420(buf, start, end, eop_lut)) {
      fprintf(stderr, "ERROR. Failed collecting EOP data\n");
      return 1;
    }
    // transform the tabulated values (time) to the TT timescale
    eop_lut.utc2tt();
    // regularize ΔUT1 and LOD values (to ΔUT1-R and LOD-R)
    eop_lut.regularize(false);
  }

  // Gravity
  // -------------------------------------------------------------------------
  // parse degree, order and the requested gravity model into a
  // HarmonicCoeffs instance. Note that to compute potential we will need
  // Lagrange polynomials (later on)
  int degree, order;
  {
    error = 0;
    error =
        dso::get_yaml_value_depth2<int>(config, "gravity", "degree", degree);
    error += dso::get_yaml_value_depth2<int>(config, "gravity", "order", order);
  }
  dso::StokesCoeffs harmonics(degree);
  {
    if (!error)
      error = dso::get_yaml_value_depth2(config, "gravity", "model", buf);
    // parse the un-normalized harmonic ceofficients from model (gfc format)
    if (!error)
      error = dso::parse_gravity_model(buf, degree, order, rnx.ref_datetime(),
                                       harmonics);
    if (error) {
      fprintf(stderr, "ERROR Failed handling gravity field model!\n");
      return 1;
    }
  }

  // Station/Beacon coordinates
  // -------------------------------------------------------------------------
  // Get beacon coordinates from sinex file and extrapolate to RINEX ref. time
  // Result coordinates per beacon are stored in the beaconCrdVec vector.
  // Note that these posotion vectors are w.r.t the beacon/antenna reference
  // point. When in actual processing, this has to be changed, if we are
  // considering iono-free analysis
  std::vector<dso::BeaconCoordinates> beaconCrdVec;
  beaconCrdVec.reserve(rnx.stations().size());
  {
    if (dso::get_yaml_value_depth2(config, "reference-frame",
                                   "station-coordinates", buf)) {
      fprintf(stderr,
              "ERROR. Failed parsing reference-frame/station-coordinates file "
              "from YAML %s\n",
              argv[1]);
      return 1;
    }
    // coordinates of beacons (on antenna RP), ECEF/pdop
    if (extrapolate_sinex_coordinates(buf, rnx.stations(), rnx.ref_datetime(),
                                      beaconCrdVec, false, true)) {
      fprintf(stderr,
              "ERROR. Failed extracting/extrapolating beacon coordinates\n");
      return 1;
    }
    printf("[DEBUG] Number of beacons in RINEX file %d, extracted coordinates "
           "for %d\n",
           (int)rnx.stations().size(), (int)beaconCrdVec.size());
  }

  /*
   * Troposphere
   * -------------------------------------------------------------------------
   * read-in the grid file. That is all for now
   */
  { dso::get_yaml_value_depth3(config, "troposphere", "vmf3", "grid", buf); }
  auto site_names_vec = beaconcrs2cchar_vec(beaconCrdVec);
  dso::SiteVMF3Feed feed(buf, site_names_vec);

  // Ocean Tides Deformation
  dso::get_yaml_value_depth2(config, "ocean-tides", "blq", buf);
  std::vector<dso::BlqSiteInfo> blqInfoVec;
  dso::Hardisp ocdeform;
  {
    // a vector containing all 4-char site names, to extract BLQ info for
    char *namepool = new char[rnx.stations().size() * 5];
    std::memset(namepool, '\0', rnx.stations().size() * 5);
    int i = 0;
    for (const auto &s : rnx.stations()) {
      std::memcpy(namepool + i * 5, s.m_station_id, sizeof(char) * 4);
      ++i;
    }
    std::vector<const char *> sites;
    sites.reserve(rnx.stations().size());
    for (int j = 0; j < (int)rnx.stations().size(); j++)
      sites.push_back(namepool + j * 5);
    if (dso::parse_blq(buf, blqInfoVec, &sites)) {
      fprintf(stderr, "Failed reading BLQ file %s\n", buf);
      return 1;
    }
    delete[] namepool;
  }

  // Ocean Tides Geopotential
  int oc_degree, oc_order;
  std::vector<dso::DoodsonOceanTideConstituent> vdds;
  dso::get_yaml_value_depth2(config, "ocean-tides", "harmonics", buf);
  error = dso::get_yaml_value_depth2<int>(config, "ocean-tides", "degree",
                                          oc_degree);
  error +=
      dso::get_yaml_value_depth2<int>(config, "ocean-tides", "order", oc_order);
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
  // An ocean tide instance
  dso::OceanTide octide(vdds, harmonics.GM(), harmonics.Re(), oc_degree,
                        oc_order);

  // Setup Integration Parameters for Orbit Integration
  // -------------------------------------------------------------------------
  // We will need the pck (SPICE) kernel for gravitational parameters of Sun
  // and Moon
  char pck_kernel[256];
  if (dso::get_yaml_value_depth2(config, "naif-kernels", "pck", pck_kernel)) {
    fprintf(stderr, "ERROR Failed locating NAIF pck kernel\n");
    return 1;
  }
  char dtm2020_data[256];
  if (dso::get_yaml_value_depth3(config, "force-model", "atmospheric-drag", "dtm2020-data", dtm2020_data)) {
    fprintf(stderr, "ERROR Failed locating dtm2020-data\n");
    return 1;
  }
  dso::IntegrationParameters IntegrationParams(degree, order, eop_lut,
                                               harmonics, pck_kernel, dtm2020_data);

  // Orbit Integrator
  // -------------------------------------------------------------------------
  // Setup an integrator, to extrapolate orbit with:
  // 1. Relative accuracy 1e-12
  // 2. Absolute accuracy 1e-12
  // 3. Num of Equations: 6 for state and 6*6 for variational equations
  //dso::SGOde Integrator(dso::VariationalEquations_thread, (6 + 6 * n), 1e-12,
  //                      1e-12, &IntegrationParams);
  dso::SGOde Integrator(dso::VariationalEquations_ta, (6 + 6 * n), 1e-12,
                        1e-12, &IntegrationParams);

  // Default observation sigma for a range-rate observable at zenith
  double sigma_obs;
  error = dso::get_yaml_value_depth2<double>(config, "filtering",
                                             "observation-sigma", sigma_obs);
  assert(sigma_obs > 0e0 && sigma_obs < 1e2 && (!error));

  // get the (RINEX) indexes for the observables we want
  int l1i, l2i, fi, w1i, w2i;
  {
    if (get_rinex_indexes(rnx, l1i, l2i, fi, w1i, w2i))
      return 1;
  }

  // Important !! OC-TIDES
  IntegrationParams.octide = &octide;
  // Setup Solid Earth Tide
  dso::SolidEarthTide setide(harmonics.GM(), harmonics.Re(),
                             IntegrationParams.GMMon * 1e9,
                             IntegrationParams.GMSun * 1e9);
  IntegrationParams.setide = &setide;

  OrbitIntegrator svState(dso::TwoPartDate(rnx.time_of_first_obs()), eop_lut);
  
  // On-board receiver eccentricity, in the satellite-fixed frame
  // -------------------------------------------------------------------------
  if (INCLUDE_ATTITUDE) {
    if (dso::get_yaml_value_depth2(config, "attitude", "mass-cog", buf)) {
      fprintf(stderr, "ERROR Failed locating Mass and CoG information file\n");
      return 1;
    }
    Eigen::Matrix<double, 3, 1> sat_cog, l3_pco;
    double sat_mass;
    // get satellite CoG coordinates in the satellite-fixed RF (with
    // corrections)
    assert(!dso::SatelliteInfo<dso::SATELLITE::Jason3>::mass_cog(
        rnx.ref_datetime(), buf, sat_mass, sat_cog));
    // get satellite ARP coordinates in the satellite-fixed RF
    Eigen::Matrix<double, 3, 1> l1_pco, l2_pco;
    dso::SatelliteInfo<dso::SATELLITE::Jason3>::pco(l1_pco, l2_pco);
    // compute the iono-free phase center
    l3_pco = l1_pco + (l2_pco - l1_pco) / (dso::GAMMA_FACTOR - 1e0);
    printf(
        "[DEBUG] Satellite CoG coordinates [%+.4f %+.4f %+.4f]m in SV-frame\n",
        sat_cog(0), sat_cog(1), sat_cog(2));
    printf("[DEBUG] AV-Antenna RP coordinates [%+.4f %+.4f %+.4f]m (Iono-Free) "
           "in SV-frame\n",
           l3_pco(0), l3_pco(1), l3_pco(2));
    // Attitude Information
    dso::get_yaml_value_depth2(config, "attitude", "body-quaternion", buf);
    /* set attitude in Integration parametrs and */
    IntegrationParams.set_sv_frame(sat_cog, l3_pco, buf, sat_mass);
    svState.set_attitude(IntegrationParams);
  }
  
  // Data feed for the DTM2020 model
  // -------------------------------------------------------------------------
  if (INCLUDE_ATM_DRAG) {
    dso::get_yaml_value_depth3(config, "force-model", "atmospheric-drag",
                               "atmo-data-csv", buf);
    const dso::TwoPartDate t0(rnx.time_of_first_obs());
    dso::SolarActivityData flux;
    if (flux.feed(buf, rnx.time_of_first_obs(), 3, 3)) {
      fprintf(stderr, "[ERROR] Failed to read required solar flux data from file %s\n", buf);
      return 1;
    }
    IntegrationParams.Dtm20.set_flux_data(flux);
  }

  /*
   * Setup the Kalman filter
   * Satellite state vector : # 6
   * m params including     : Cd, Csrp
   * Δfe/feN per beacon     : # NumBeacons
   */
  const int NumBeacons = beaconCrdVec.size();
  const int NumParams = 6 + m + NumBeacons;
  double apriori_sigma_freqbias, apriori_sigma_orbsp3_pos,
      apriori_sigma_orbsp3_vel, apriori_sigma_Cd, apriori_sigma_Cr;
  Kalman filter(NumParams);
  { /* a-priori P matrix */
    error = dso::get_yaml_value_depth2<double>(
        config, "filtering", "deffen-sigma", apriori_sigma_freqbias);
    error = dso::get_yaml_value_depth2<double>(
        config, "filtering", "orbsp3-sigma-pos", apriori_sigma_orbsp3_pos);
    error = dso::get_yaml_value_depth2<double>(
        config, "filtering", "orbsp3-sigma-vel", apriori_sigma_orbsp3_vel);
    error = dso::get_yaml_value_depth2<double>(config, "filtering", "cd-sigma",
                                               apriori_sigma_Cd);
    error = dso::get_yaml_value_depth2<double>(config, "filtering", "cr-sigma",
                                               apriori_sigma_Cr);
    /* Initialize to Unit matrix */
    filter.P = Eigen::MatrixXd::Identity(NumParams,NumParams);
    /* Covariance for state/position */
    filter.P.block<3, 3>(0, 0) *=
        (apriori_sigma_orbsp3_pos * apriori_sigma_orbsp3_pos);
    /* Covariance for state/velocity */
    filter.P.block<3, 3>(3, 3) *=
        (apriori_sigma_orbsp3_vel * apriori_sigma_orbsp3_vel);
    /* Covariance for Δfe / feN */
    filter.P.block(6 + m, 6 + m, NumBeacons, NumBeacons) *=
        (apriori_sigma_freqbias * apriori_sigma_freqbias);
    /* drag coefficient */
    if (m >= 1) {
      filter.P(6, 6) *= (apriori_sigma_Cd * apriori_sigma_Cd);
      if (m > 1)
        /* radiation coefficient */
        filter.P(7, 7) *= (apriori_sigma_Cr * apriori_sigma_Cr);
    }
  }
  double apriori_Cd, apriori_Cr;
  { /* a-priori estimates */
    if (m >= 1) {
      dso::get_yaml_value_depth3<double>(config, "force-model",
                                         "atmospheric-drag", "Cd-apriori", apriori_Cd);
      dso::get_yaml_value_depth3<double>(config, "force-model", "srp",
                                         "Cr-apriori", apriori_Cr);
      filter.x(6) = apriori_Cd;
      if (m > 1) {
        filter.x(7) = apriori_Cr;
      }
    }
    for (int i=6+m; i<NumParams; i++) filter.x(i) = 1e-3;
  }

  // Start RINEX data-block iteration
  // -------------------------------------------------------------------------
  // get an iterator to the RINEXs data blocks
  dso::RinexDataBlockIterator it(&rnx);

  // Some variables ...
  [[maybe_unused]] const double J2 = harmonics.J2();
  [[maybe_unused]] const double GM = harmonics.GM();
  [[maybe_unused]] const double Re = harmonics.Re();
  // form a rotation instance (ITRF-to-GCRF) for current epoch
  dso::Itrs2Gcrs Rot(rnx.time_of_first_obs(), &eop_lut);
  // a buffer to write datetime strings to ...
  char dtbuf[64];
  // counters ...
  [[maybe_unused]] unsigned flaged_obs =
      0; // number of observations with 'bad' flags
  [[maybe_unused]] unsigned num_obs =
      0; // observation count, regardless if usable or not
  [[maybe_unused]] unsigned num_blocks = 0; // RINEX block count
  /* number of observations with low elevation (rejected) */
  [[maybe_unused]] unsigned flagged_low_obs = 0;
  /* number of observations flagged as outliers */
  [[maybe_unused]] unsigned flagged_res_obs = 0;
  /* number of observations flagged for SAA */
  [[maybe_unused]] unsigned flagged_saa_obs = 0;
  /* number of observations used for filtering */
  [[maybe_unused]] unsigned filter_obs = 0;
  /* currently in SAA */
  [[maybe_unused]] int inside_saa = 0;
  // error flag
  error = 0;

  /* Epoch when the last revolution was concluded */
  dso::TwoPartDate last_revolution_at(dso::datetime<dso::nanoseconds>::max());

  /* store last beacon observation */
  std::vector<RcSv> vLastObs;
  vLastObs.reserve(rnx.stations().size());

  // for every new data block in the RINEX file (aka every epoch) ...
  while (!(error = it.next())) {
    // current proper time (aka τ)
    [[maybe_unused]] const auto tobs_proper =
        dso::TwoPartDate(it.proper_time());
    [[maybe_unused]] const auto tobs_proper_dt = it.proper_time();

    // get current observation epoch (tobs) from proper time to TAI
    [[maybe_unused]] const auto tobs_tai_dt = it.corrected_l1_epoch();
    [[maybe_unused]] const auto tobs_tai = dso::TwoPartDate(tobs_tai_dt);

    // get current observation epoch in UTC
    [[maybe_unused]] const auto tobs_utc = tobs_tai.tai2utc();

    /* have the current date (TAI) ready to print in dtbuf */
    dso::strftime_ymd_hmfs(tobs_tai_dt, dtbuf);

    // integrate orbit to here (TAI)
    if (!num_blocks) {
      // for first iteration, get reference state from sp3, for an epoch as
      // close as possible
      if (sp3_iterator.goto_epoch(tobs_tai_dt)) {
        fprintf(stderr, "ERROR Failed to get reference position from SP3\n");
        return 1;
      }
      Eigen::Matrix<double, 6, 1> itrf;
      itrf(0) = sp3_iterator.data_block().state[0] * 1e3;
      itrf(1) = sp3_iterator.data_block().state[1] * 1e3;
      itrf(2) = sp3_iterator.data_block().state[2] * 1e3;
      itrf(3) = sp3_iterator.data_block().state[4] * 1e-1;
      itrf(4) = sp3_iterator.data_block().state[5] * 1e-1;
      itrf(5) = sp3_iterator.data_block().state[6] * 1e-1;
      /* state in ITRF, w.r.t CoM */
      svState.set_state_from_itrf(dso::TwoPartDate(sp3_iterator.current_time()),
                                  itrf);
      if (svState.integrate(tobs_tai, Integrator)) {
        fprintf(stderr, "ERROR. Failed to integrate orbit!\n");
        return 1;
      }
      /* assign time of initialization; considered time of first revolution */
      last_revolution_at = tobs_tai;
    } else {
      if (svState.integrate(tobs_tai, Integrator)) {
        fprintf(stderr, "ERROR. Failed to integrate orbit!\n");
        return 1;
      }
    }

    /* an ITRF-to-GCRF Rotation for general use */
    dso::Itrs2Gcrs R(tobs_tai.tai2tt(), &eop_lut);

    /* filter time-update */
    {
      Eigen::MatrixXd F = Eigen::MatrixXd::Identity(NumParams, NumParams);
      F.block<6, 6>(0, 0) = svState.stateTransitionMatrix();
      filter.time_update(tobs_tai, F);
      filter.estimates().block<6, 1>(0, 0) = svState.gcrf_state_cm();
      if (m>=1) IntegrationParams.drag_ceofficient() = filter.x(6);
      if (m>1) IntegrationParams.srp_ceofficient() = filter.x(7);
    }

    /* check to see if we have a new orbit revolution */
    if (tobs_tai.diff<dso::DateTimeDifferenceType::FractionalSeconds>(
            last_revolution_at) > 112.42 * 60e0) {
      printf("[REV] New SV revolution around the Earth at %s "
             "re-initialize Cd\n",
             dtbuf);
      //if (m >= 1) {
      //  // filter.x(6) = apriori_Cd;
      //  filter.P.row(6).setZero();
      //  filter.P.col(6).setZero();
      //  filter.P(6, 6) = (apriori_sigma_Cd * apriori_sigma_Cd);
      //}
      //if (m >= 1)
      //  IntegrationParams.drag_ceofficient() = filter.x(6);
      last_revolution_at = tobs_tai;
    }

    /* block-specific counters */
    int num_cblock_obs = 0, flagged_low_obs_block = 0,
        flagged_saa_obs_block = 0, flagged_res_obs_block = 0,
        filter_obs_block = 0;
    inside_saa = 0;

    /* iterate through beacons in current block (epoch) */
    auto beaconIt = it.cblock.begin();
    while (beaconIt != it.cblock.end()) {
      /* augment counters */
      ++num_obs;
      ++num_cblock_obs;
      /* Beacon's 4-char id */
      const char *b4id = rnx.beacon_internal_id2id(beaconIt->id());
      /* check if beacon is within SSA */
      if (!site_affected_by_ssa(b4id)) {
      /* check flags */
      if (check_obs_flags(*beaconIt, l1i, l2i, w1i, w2i)) {
        ++flaged_obs;
        /* mark discontinuity */
        if (auto vit = findLastObs(b4id, vLastObs); vit != vLastObs.end())
          vit->mark_restart();
      } else { /* observation falgs ok, keep on */
        /* index for referencing beacon parameters in filter */
        int beaconFilterIndex;
        /* coordinates of beacon, ITRF */
        Eigen::Matrix<double, 3, 1> bcrd;
        if (auto vit = findItrfCrd(b4id, beaconCrdVec);
            vit != beaconCrdVec.end()) {
          bcrd = getItrfCrd(vit);
          beaconFilterIndex = std::distance(beaconCrdVec.cbegin(), vit) + 6 + m;
        /* Rc-Sv azimouth [rad], elevation [rad] and range [m] */
        double az, el;
        {
          Eigen::Matrix<double, 3, 1> r_enu =
              dso::car2top<dso::ellipsoid::grs80>(bcrd,
                                                  svState.itrf_position_arp());
          dso::top2dae(r_enu, az, el);
        }
        /* only procced if elevation angle is large enough */
        if (dso::rad2deg(el) > EleCutOff) {
          /* get beacon nominal frequency, f_eN at 2GHz in [Hz] */
          double feN;
          {
            int k; // shift factor
            if (rnx.beacon_shift_factor(beaconIt->id(), k)) {
              fprintf(stderr,
                      "[ERROR] Failed to find shift factor for beacon %.4s\n",
                      b4id);
              return 1;
            }
            double fe2N; // [Hz]
            dso::beacon_nominal_frequency(k, feN, fe2N);
          }
          /* Iono correction in [cycles] */
          const double Dion = iono_l2_correction(*beaconIt, l1i, l2i);
          /* Tropospheric correction (only proceed if ok) */
          TropoDetails Dtro;
          if (get_tropo_vmf(b4id, tobs_tai_dt, bcrd, dso::DPI / 2e0 - el, feed,
                            Dtro)) {
            fprintf(stderr,
                    "[WRNNG] Failed to get tropo ceorrection for %.4s at %s; "
                    "observation skipped\n",
                    b4id, dtbuf);
            // return 1;
          } else {
            /* Hold current measurement details in a struct */
            RcSv cObs(b4id, tobs_tai, beaconIt->m_values[l1i].m_value, Dion,
                      Dtro, svState.gcrf_position_arp(), bcrd, eop_lut);
            /* Find emission time (correction of aberration */
            {
              const dso::TwoPartDate tai_emission =
                  correction_aberration(cObs, R);
              cObs.set_emission_tai(tai_emission, eop_lut);
            }
            /* set proper time */
            cObs.set_proper_time(tobs_proper);
            /* depending on if we already have a records for the beacon */
            const auto pObs = std::find_if(
                vLastObs.begin(), vLastObs.end(), [=](const RcSv &v) {
                  return !(std::strncmp(b4id, v.bc_name(), 4));
                });
            if (pObs == vLastObs.end()) {
              // cObs.mark_restart();
              vLastObs.emplace_back(cObs);
            } else {
              if (/*pObs->restart() ||*/
                  tobs_tai.diff<dso::DateTimeDifferenceType::FractionalSeconds>(
                      pObs->time_of_reception()) > RESTART_AFTER_SEC) {
                /* restart frequency bias for new pass if needed */
                if (tobs_tai
                        .diff<dso::DateTimeDifferenceType::FractionalSeconds>(
                            pObs->time_of_reception()) >
                    NEW_PASS_AFTER_MIN * 60e0) {
                  filter.reset_frequency_bias(
                      beaconFilterIndex, filter.estimates()(beaconFilterIndex),
                      apriori_sigma_freqbias);
                  printf("[DEBUG] Resetting frequency bias for site %.4s, new "
                         "pass at %s\n",
                         b4id, dtbuf);
                
                  cObs.mark_restart();
                  *pObs = cObs;
                }
              } else {
                /* observation equation */
                double Vobs, Vtheo;
                /* true proper frequency of the receiver f_rT [Hz] */
                const double frT =
                    dso::DORIS_FREQ1_MHZ * 1e6 *
                    (1e0 + rfo_fit.value_at(tobs_proper_dt) * 1e-11);
                const double DfeFen = filter.parameter(beaconFilterIndex);
                /* emitter potential */
                const double ePot =
                    emitter_potential(bcrd, R.omega_earth(), 3.986004418e14);
                /* partials */
                double d_DfedeN;
                Eigen::Matrix<double, 3, 1> dObsdr, dObsdv;
                observation_equation(*pObs, cObs, feN, frT, DfeFen, ePot, Vobs,
                                     Vtheo, d_DfedeN, dObsdr, dObsdv);
                /* filter observation update */
                {
                  Eigen::VectorXd H = Eigen::VectorXd::Zero(NumParams);
                  H.block<3, 1>(0, 0) = dObsdr;
                  H.block<3, 1>(3, 0) = dObsdv;
                  if (m>=1) H(6) = 1e0;
                  if (m>1) H(7) = 1e0;
                  H(beaconFilterIndex) = d_DfedeN;
                  double var_prediction;
                  double res_prediction = filter.prediction_residual(
                      Vobs, -Vtheo, sigma_obs, H, var_prediction);
                  /* outlier check */
                  if (std::abs(res_prediction) > 3 * std::sqrt(var_prediction)) {
                    ++flagged_res_obs;
                    ++flagged_res_obs_block;
                  } else if (cObs.restart()) {
                    cObs.mark_restart(0);
                    filter.observation_no_update(Vobs, -Vtheo, sigma_obs / std::sin(el) / std::sin(el), H);
                  } else {
                    filter.observation_update(Vobs, -Vtheo,
                                              sigma_obs / std::sin(el), H);
                    /* update drag coefficient value */
                    if (m>=1) IntegrationParams.drag_ceofficient() = filter.x(6);
                    if (m>1) IntegrationParams.srp_ceofficient() = filter.x(7);
                    /* debug print */
                    printf("[RES] %s site:%.4s res:%+.6f res_prediction:%+.6f "
                           "res_var:%.6f vobs:%+.3f vtheo:%+.3f df:%.3e "
                           "df_var:%.9e cd:%.6f cd_var:%.9e cr:%.6f "
                           "cr_var:%.9e el:%.3f\n",
                           dtbuf, b4id, Vobs + Vtheo, res_prediction,
                           var_prediction, Vobs, Vtheo, DfeFen,
                           filter.P(beaconFilterIndex, beaconFilterIndex),
                           IntegrationParams.drag_ceofficient(), filter.P(6, 6),
                           IntegrationParams.srp_ceofficient(), filter.P(7, 7),
                           dso::rad2deg(el));
                    ++filter_obs;
                    ++filter_obs_block;
                  }
                }
                /* push back */
                *pObs = cObs;
              }
            } /* tropospheric correction found/applied */
          }   /* end handling current observation, cObs */
        }   /* end of observation with acceptable elevation */ 
        else if (dso::rad2deg(el) < 0) {
          fprintf(
              stderr,
              "[WRNNG] Unexpected elevation angle encountered for %.4s at %s\n",
              b4id, dtbuf);
        } else {
          ++flagged_low_obs;
          ++flagged_low_obs_block;
        }
        } else {
          fprintf(stderr,
                  "[ERROR] Failed to find coordinates for station %4s\n", b4id);
        }
      } /* end of non-flaged observation processing */
    } else {
      // fprintf(stderr, "[WRNNG] Site %.4s within SSA; ignornig observation\n", b4id);
      ++flagged_saa_obs;
      ++flagged_saa_obs_block;
      inside_saa = 1;
    } /* site affected by SSA */
      ++beaconIt;
    } /* end of beacons in current block */
      /* print block statistics */
    printf("[BLC] %s obs:%d rele:%d rsaa:%d rres:%d ufilter:%d\n", dtbuf,
           num_cblock_obs, flagged_low_obs_block, flagged_saa_obs_block,
           flagged_res_obs_block, filter_obs_block);

    { // print state
      printf("[ORB] %s x:%+.6f y:%+.6f z:%+.6f vx:%+.9f vy:%+.9f vz:%+.9f "
             "x_var:%+.6e y_var:%+.6e z_var:%+.6e vx_var:%+.9e vy_var:%+.9e "
             "vz_var:%+.e\n",
             dtbuf, svState.itrf_state_cm()(0), svState.itrf_state_cm()(1),
             svState.itrf_state_cm()(2), svState.itrf_state_cm()(3),
             svState.itrf_state_cm()(4), svState.itrf_state_cm()(5),
             filter.P(0, 0), filter.P(1, 1), filter.P(2, 2), filter.P(3, 3),
             filter.P(4, 4), filter.P(5, 5));
      printf("[ECI] %s x:%+.9f y:%+.9f z:%+.9f vx:%+.12e vy:%+.12e vz:%+.12e\n",
             dtbuf, svState.gcrf_state_cm()(0), svState.gcrf_state_cm()(1),
             svState.gcrf_state_cm()(2), svState.gcrf_state_cm()(3),
             svState.gcrf_state_cm()(4), svState.gcrf_state_cm()(5));
    }

    ++num_blocks; /* augment data block counter */

    if (tobs_tai.diff<dso::DateTimeDifferenceType::FractionalDays>(
            dso::TwoPartDate(rnx.time_of_first_obs())) > MAX_HOURS / 24e0)
      break;
  } /* data blocks ended, no more data in RINEX */

  printf("#[DEBUG] Number of data block read: %u\n", num_blocks);
  return 0;
}
