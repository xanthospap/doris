#include "doris_rinex.hpp"
#include <algorithm>

// Reminder: in the RINEX files, the F value (aka relative frequency offset for
// the receiver's oscillator) is given in units: 1e-11

class LinearKalman {
public:
  LinearKalman(const double *state0, const double *state_stddev, double snoise,
               const dso::datetime<dso::nanoseconds> &reft) noexcept {
    state[0] = state0[0];
    state[1] = state0[1];
    P[0] = state_stddev[0] * state_stddev[0]; // P(1,1)
    P[1] = state_stddev[0] * state_stddev[1]; // P(1,2) = P(2,1)
    P[2] = state_stddev[1] * state_stddev[1]; // P(2,2)
    R = snoise * snoise;
    tref = reft;
  }

  double model_value(const dso::datetime<dso::nanoseconds> &t) const noexcept {
    double dt = sdt(t);
    return state[0] + dt * state[1];
  }

  double sdt(const dso::datetime<dso::nanoseconds> &t) const {
    auto nsec = t.delta_sec(tref);
    // return nsec.to_fractional_seconds();
    return nsec.fractional_days();
  }

  void update(double z, const dso::datetime<dso::nanoseconds> &t) noexcept {
    double dt = sdt(t);

    /* model y = a + b * t
     * Notes:
     * measurement z(k) := a + b *t => H = [1 t]
     * model: y(k+1) = y(k) (a and b are constant), aka F=I(2x2)
     */
    F = np.eye(2, 2)
    H = np.mat([ 1e0, dt ])
    newx = F *self.x 
    newP = F * self.P * F.transpose()
    K = newP * H.transpose() / (H * newP * H.transpose() + self.R) 
    self.x = newx + K * (z - H * newx) 
    self.P = (np.eye(2) - K * H) *newP + self.Q

    // measurement pre-fit residual
    double y = z - state[0];
    // pre-fit residual covariance
    double s = P[0] + R;
    assert(s != 0e0);
    // kalman gain
    K[0] = P[0] / s;
    K[1] = P[2] / s;
    // updated state estimate
    state[0] += K[0] * y;
    state[1] += K[1] * y;
    // update estimate covariance
    double p0 = (1e0 - K[0]) * P[0];
    double p1 = (1e0 - K[0]) * P[1];
    double p2 = P[2];
    P[0] = p0;
    P[1] = p1;
    P[2] = p2;
    // predicted state estimate
    state[0] = state[0] + state[1] * dt;
    state[1] = state[1];
    // predicted estimate covariance
    p0 = P[0] + dt * P[1] + dt * (P[1] + dt * P[2]);
    p1 = P[1] + P[2] * dt;
    p2 = P[2];
    P[0] = p0; //+ q[0];
    P[1] = p1; //+ q[1];
    P[2] = p2; //+ q[2];
    return;
  }

  double state[2];
  double R;
  double K[2];
  double P[3];
  dso::datetime<dso::nanoseconds> tref;
};

int get_starting_rfo(const char *rnx_fn, double &rfo) noexcept {
  using namespace ids;

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

  rnx.stream().getline(line, DorisObsRinex::MAX_RECORD_CHARS);
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
    fprintf(stderr,
            "[ERROR] Failed to resolve data block! error=%d (traceback: %s)\n",
            status, __func__);
    return 1;
  }

  // for a single block, all F measurements/values should be the same!
  rfo = obsvec[0].m_values[findex].m_value;

  // return, all ok
  return 0;
}

int get_rinex_rfo(const char *rnx_fn, long &rfo_collected,
                  LinearKalman *filter) {
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

  double rfo = 0e0;
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
    rfo = obsvec[0].m_values[findex].m_value;
    for (const auto &bobs : obsvec) {
      if (bobs.m_values[findex].m_value != rfo) {
        fprintf(stderr,
                "[ERROR] F value differs within the same data block! "
                "(traceback: %s)\n",
                __func__);
        return 1;
      }
    }

    // update via kalman using this new F value
    filter->update(rfo, hdr.m_epoch);
    printf("%.8f %.8f %.8f +/- %.10f %.8f +/- %.10f\n", hdr.m_epoch.as_mjd(),
           rfo, filter->state[0], std::sqrt(filter->P[0]), filter->state[1],
           std::sqrt(filter->P[2]));
    ++idx;
  }

  rfo_collected = idx;

  if (rnx.stream().eof()) {
    rnx.stream().clear();
    return 0;
  }

  return 1;
}

int ids::fit_relative_frequency_offset(char **rinex_fns, int num_rinex,
                                       double sigma_x,
                                       double sigma_vx,
                                       double sigma_z) noexcept {

  long rfo_collected = 0;
  int error = 0;

  // let's get an initial estimate for RFO using the first measurement in the
  // first RINEX file
  double rfo_0 = 0e0;
  try {
    if (get_starting_rfo(rinex_fns[0], rfo_0)) {
      fprintf(stderr,
              "[ERROR] Failed to read F value from RINEX %s (traceback: %s)\n",
              rinex_fns[0], __func__);
      return 1;
    }
  } catch (std::exception &e) {
    fprintf(
        stderr,
        "[ERROR] Exception caught while reading first F value; what string: "
        "%s (traceback: %s)\n",
        e.what(), __func__);
    return 1;
  }

  double x[] = {rfo_0, 1e-1};
  double xstd[] = {sigma_x, sigma_vx};
  double snoise = sigma_z;
  LinearKalman kalman(x, xstd, snoise,
                      dso::datetime<dso::nanoseconds>{dso::year(2021),
                                                      dso::month(1),
                                                      dso::day_of_month(1)});

  // for every rinex in the list ...
  for (int i = 0; i < num_rinex; i++) {
    try {
      error = get_rinex_rfo(rinex_fns[i], rfo_collected, &kalman);
      if (error)
        return error;
    } catch (std::exception &e) {
      fprintf(stderr,
              "[ERROR] Exception caught while fitting F values; what string: "
              "%s (traceback: %s)\n",
              e.what(), __func__);
      fprintf(stderr,
              "[ERROR] Failed to collectd/update F values for RINEX %s "
              "(traceback: %s)\n",
              rinex_fns[i], __func__);
      return 1;
    }
  }

  // just for debuging ... reconstruct the model
  auto t1 = dso::datetime<dso::nanoseconds>{dso::year(2021), dso::month(1),
                                            dso::day_of_month(3)};
  auto t2 = dso::datetime<dso::nanoseconds>{dso::year(2021), dso::month(1),
                                            dso::day_of_month(4)};
  while (t1 < t2) {
    printf("Fit: %.8f %.8f\n", t1.as_mjd(), kalman.model_value(t1));
    t1.add_seconds(dso::seconds(3));
  }

  return 0;
}