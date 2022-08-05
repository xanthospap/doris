#include "doris_rinex.hpp"
#include <algorithm>

// Reminder: in the RINEX files, the F value (aka relative frequency offset for
// the receiver's oscillator) is given in units: 1e-11

struct LinearKalman {
  LinearKalman(const double *x0, const double *x0_stddev, double obs_stddev,
               const dso::datetime<dso::nanoseconds> &reft) noexcept {
    x[0] = x0[0];
    x[1] = x0[1];
    P[0] = x0_stddev[0] * x0_stddev[0]; // P(1,1)
    P[1] = 0e0; // x0_stddev[0] * x0_stddev[1]; // P(1,2) = P(2,1)
    P[2] = x0_stddev[1] * x0_stddev[1]; // P(2,2)
    R = obs_stddev * obs_stddev;
    tref = reft;
  }

  double value_at(const dso::datetime<dso::nanoseconds> &t) const noexcept {
    double dt = t2d(t);
    return x[0] + dt * x[1];
  }

  double t2d(const dso::datetime<dso::nanoseconds> &t) const {
    auto nsec = t.delta_sec(tref);
    return nsec.fractional_days();
  }

  void update(double z, const dso::datetime<dso::nanoseconds> &t) noexcept {
    double dt = t2d(t);

    double y = z - value_at(t);
    double p0p1dt = P[0] + P[1]*dt;
    double p1p2dt = P[1] + P[2]*dt;
    double s = p0p1dt + dt*p1p2dt + R;
    double k0 = p0p1dt / s;
    double k1 = p1p2dt / s;
    x[0] += k0*y;
    x[1] += k1*y;
    double p0 = (1e0-k0)*P[0] - P[1]*k0*dt;
    double p1 = ( ((1e0-k0)*P[1] - P[2]*k0*dt) +
      (-k1*P[0] + (1e0-k1*dt)*P[1]) ) / 2e0;
    double p2 = -k1*P[1] + (1e0-k1*dt) * P[2];
    P[0] = p0;
    P[1] = p1;
    P[2] = p2;

    /* predict
    x[0] += dt * x[1];
    x[1] = x[1];

    double p0p1dt = P[0] + P[1]*dt;
    double p1p2dt = P[1] + P[2]*dt;
    P[0] = p0p1dt + dt*p1p2dt + Q[0];
    P[1] = p1p2dt;
    P[2] = P[2] + Q[1];

    // update
    double y = z - value_at(t);
    double s = P[0] + R;
    double k0 = P[0] / s;
    double k1 = P[1] / s;
    x[0] += k0 * y;
    x[1] += k1 * y;

    double p0 = P[0] * (1e0-k0);
    double p1 = (P[1] * (1e0-k0) + P[1] - k1*P[0]) / 2e0;
    printf(">> diff =%.8f\n", P[1] * (1e0 - k0) - P[1] + k1 * P[0]);
    double p2 = P[2] - k1*P[1];

    P[0] = p0;
    P[1] = p1;
    P[2] = p2;
    */

    return;
  }

  double x[2];
  double K[2];
  double P[3];
  double R;
  double Q[2] = {0e0, 0e0};
  dso::datetime<dso::nanoseconds> tref;
};

int get_starting_rfo(const char *rnx_fn, double &rfo,
                     dso::datetime<dso::nanoseconds> &reft) noexcept {
  using namespace dso;

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

  // get & assign the header time
  reft = hdr.m_epoch;

  // return, all ok
  return 0;
}

int get_rinex_rfo(const char *rnx_fn, long &rfo_collected,
                  LinearKalman *filter) {
  using namespace dso;

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
    if (int error = rnx.read_data_block(hdr, obsvec); error) {
      fprintf(
          stderr,
          "[ERROR] Failed to resolve data block! error=%d (traceback: %s)\n",
          error, __func__);
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
           rfo, filter->x[0], std::sqrt(filter->P[0]), filter->x[1],
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

int dso::fit_relative_frequency_offset(char **rinex_fns, int num_rinex,
                                       double sigma_x, double sigma_vx,
                                       double sigma_z) noexcept {

  long rfo_collected = 0;
  int error = 0;

  // let's get an initial estimate for RFO using the first measurement in the
  // first RINEX file; also get reference time for kalman
  dso::datetime<dso::nanoseconds> reft;
  double rfo_0 = 0e0;
  try {
    if (get_starting_rfo(rinex_fns[0], rfo_0, reft)) {
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

  // initialize the (linear) kalman filter
  double x[2] = {rfo_0, 3e0};
  double xstd[2] = {sigma_x, sigma_vx};
  LinearKalman kalman(x, xstd, sigma_z, reft);

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
                                            dso::day_of_month(1)};
  auto t2 = dso::datetime<dso::nanoseconds>{dso::year(2021), dso::month(1),
                                            dso::day_of_month(8)};
  while (t1 < t2) {
    printf("Fit: %.8f %.8f\n", t1.as_mjd(), kalman.value_at(t1));
    t1.add_seconds(dso::seconds(30));
  }
  printf("Model: %.4f + %.4f * x\n", kalman.x[0], kalman.x[1]);

  return 0;
}
