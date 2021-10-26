#include "doris_rinex.hpp"
#include <algorithm>

// Reminder: in the RINEX files, the F value (aka relative frequency offset for
// the receiver's oscillator) is given in units: 1e-11

class LinearKalman {
public:
  LinearKalman(const double *state0, const double *state_stddev, double snoise, const dso::datetime<dso::nanoseconds> &reft) noexcept {
    state[0] = state0[0];
    state[1] = state0[1];
    P[0] = state_stddev[0] * state_stddev[0];
    P[1] = state_stddev[0] * state_stddev[1];
    P[2] = state_stddev[1] * state_stddev[1];
    R = snoise * snoise;
    tref = reft;
  }

  double sdt(const dso::datetime<dso::nanoseconds> &t) const {
    auto nsec = t.delta_sec(tref);
    // return nsec.to_fractional_seconds();
    return nsec.fractional_days();
  }

  void update(double z, const dso::datetime<dso::nanoseconds> &t) noexcept {
    double dt = sdt(t);
    // printf(">> dt=%20.10f\n", dt);

    // State prediction
    double xk0 = state[0] + state[1] * dt;
    double xk1 = state[1];
    state[0] = xk0;
    state[1] = xk1;
    
    P[0] = P[0] + dt*P[1]+dt*(P[1]+dt*P[2]);
    P[1] = P[1] + dt*P[2];
    P[2] = P[2];
    
    // measurement prediction
    double s = P[0] + R;
    assert(s!=0e0);
    double zk = state[0];

    // update
    K[0] = P[0] / s;
    K[1] = P[1] / s;
    // printf(">> update dt=%.5f dz(%.5f-%.5f)=%.5f K=[%.5f, %.5f]\n", dt, z, zk, (z-zk), K[0], K[1]);
    
    xk0 = state[0] + K[0] * (z-zk);
    xk1 = state[1] + K[1] * (z-zk);
    state[0] = xk0;
    state[1] = xk1;

    // printf(">> update P: P[0]=%.5f P[1]=%.5f P[2]=%.5f, s=%.5f\n", P[0], P[1], P[2], s);
    double P0 = P[0] - K[0]*K[0]*s;
    double P1 = P[1] - K[0]*K[1]*s;
    double P2 = P[2] - K[1]*K[1]*s;
    P[0] = P0;
    P[1] = P1;
    P[2] = P2;

    return;
  }

  double state[2];
  double R;
  double K[2];
  double P[3];
  dso::datetime<dso::nanoseconds> tref;
};

int get_rinex_rfo(const char *rnx_fn, long &rfo_collected, LinearKalman *filter) {
  using namespace ids;

  rfo_collected = 0;
  DorisObsRinex rnx(rnx_fn); // may throw ....
  char line[DorisObsRinex::MAX_RECORD_CHARS];
  RinexDataRecordHeader hdr;
  std::vector<BeaconObservations> obsvec;

  // index of observable F in the list of observation codes for this RINEX
  auto obs_list = rnx.observation_codes();
  auto fit = std::find_if(
      obs_list.begin(), obs_list.end(), [](const ObservationCode &o) {
        return o.type() == ObservationType::frequency_offset;
      });
  if (fit == obs_list.end()) {
    fprintf(stderr,
            "[ERROR] Failed to find observation type F in list of RINEX "
            "observables! (traceback: %s)\n",
            __func__);
    return 1;
  }
  int findex = std::distance(obs_list.begin(), fit);

  int idx = 0;
  while (rnx.stream() &&
         rnx.stream().getline(line, DorisObsRinex::MAX_RECORD_CHARS)) {

    // resolve the data-block header
    if (int status = rnx.resolve_data_epoch(line, hdr); status) {
      fprintf(stderr,
              "[ERROR] Failed to resolve data block header! error=%d "
              "(traceback: %s)\n",
              status, __func__);
      return 1;
    }

    // get the measurements
    if (int status = rnx.read_data_block(hdr, obsvec); status) {
      fprintf(
          stderr,
          "[ERROR] Failed to resolve data block! error=%d (traceback: %s)\n",
          status, __func__);
      return 1;
    }

    // for a single block, all F measurements/values should be the same!
    double rfo = obsvec[0].m_values[findex].m_value;
    for (const auto &bobs : obsvec) {
      if (bobs.m_values[findex].m_value != rfo) {
        fprintf(stderr,
                "[ERROR] F value differs within the same data block! "
                "(traceback: %s)\n",
                __func__);
        return 1;
      }
    }

    // push back the reolved time and offest
    // rfo[offset + idx] = rfo;
    // t[offset + idx] = hdr.m_epoch;
    filter->update(rfo, hdr.m_epoch);
    printf("%.8f %.8f %.8f +/- %.10f %.8f +/- %.10f\n", hdr.m_epoch.as_mjd(), rfo, filter->state[0], std::sqrt(filter->P[0]), filter->state[1], std::sqrt(filter->P[2]));
    ++idx;
  }

  rfo_collected = idx;

  if (rnx.stream().eof()) {
    rnx.stream().clear();
    return 0;
  }

  return 1;
}

int ids::fit_relative_frequency_offset(char **rinex_fns, int num_rinex) noexcept {

  long rfo_collected = 0;
  int error = 0;
  
  double x[] = {7.5e3, 1e-3};
  double xstd[] = {5e1, 1e-1};
  double snoise = 1e1;
  LinearKalman kalman(x, xstd, snoise, dso::datetime<dso::nanoseconds>{dso::year(2021), dso::month(1), dso::day_of_month(1)});

  // for every rinex in the list ...
  for (int i = 0; i < num_rinex; i++) {
    try {
      error = get_rinex_rfo(rinex_fns[i], rfo_collected, &kalman);
      if (error) return error;
    } catch (std::exception &e) {
      fprintf(stderr, "[ERROR] Exception caught while fitting F values; what string: %s (traceback: %s)\n", e.what(), __func__);
      fprintf(stderr,
              "[ERROR] Failed to collectd/update F values for RINEX %s "
              "(traceback: %s)\n",
              rinex_fns[i], __func__);
      return 1;
    }
  }

  // just for debuging ... reconstruct the model
  auto t1 = dso::datetime<dso::nanoseconds>{dso::year(2021), dso::month(1), dso::day_of_month(3)};
  auto t2 = dso::datetime<dso::nanoseconds>{dso::year(2021), dso::month(1), dso::day_of_month(6)};
  while (t1<t2) {
    auto nsec = t1.delta_sec(kalman.tref);
    double dt = nsec.fractional_days();
    double val = kalman.state[0] + kalman.state[1] * dt;
    printf("Fit: %.8f %.8f %.8f\n", t1.as_mjd(), val, dt);
    t1.add_seconds(dso::seconds(3));
  }

  return 0;
}