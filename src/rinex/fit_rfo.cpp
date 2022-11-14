#include "doris_rinex.hpp"
//#include "datetime/dtcalendar.hpp"
#include <algorithm>
#include <numeric>

typedef dso::datetime<dso::nanoseconds> Datetime;

// Reminder: in the RINEX files, the F value (aka relative frequency offset for
// the receiver's oscillator) is given in units: 1e-11

int parse_rinex_rfos(const char *fnrnx, 
    bool use_tai,
    int every,
    Datetime &tref,
    const dso::PolynomialModel<Datetime> &fit,
    Eigen::VectorXd &b,
    Eigen::MatrixXd &A,
    unsigned &idx,
    const Datetime &start,
    const Datetime &end) noexcept {

  // polynomial order
  const int poly_order = fit.order();
  
  // Construct the instance
  dso::DorisObsRinex rnx(fnrnx);

  // Set reference date (if not set)
  if (tref == Datetime::min())
    tref = rnx.ref_datetime();
  
  // index of the F measurement (relative frequency offset)
  const int f_idx = rnx.get_observation_code_index(
      dso::ObservationCode{dso::ObservationType::frequency_offset});

  // get an iterator to the RINEXs data blocks
  dso::RinexDataBlockIterator it(&rnx);

  // buffer for F values
  std::vector<double> fs(every,0e0);

  int error, cval=0;
  dso::datetime<dso::nanoseconds> t0, tend;
  // for every new data block in the RINEX file (aka every epoch) ...
  while (!(error = it.next())) {
    // the current reference time for the L1 observation (corrected for
    // receiver clock offset). That is tl1 is approximately TAI.
    const auto tl1 = it.corrected_l1_epoch();

    if (tl1 >= start && tl1 < end) {

      // current proper time (aka tau)
      const auto tproper = it.proper_time();

      // get the first observation set; note that all observations in block
      // should have the same F value
      auto beaconobs = it.cblock.begin();

      // temporarily store value at buffer
      fs[cval % every] = beaconobs->m_values[f_idx].m_value;

      // store current date
      tend = (use_tai) ? tl1 : tproper;
      
      // check if we reached every values
      if (!((cval + 1) % every)) {
        // we did, compute data points and add to matrices
        b(idx) = std::reduce(fs.begin(), fs.end(), 0e0) / every;

        // midpoint in time, between t0 and tend
        const auto dsec = tend.delta_sec(t0);
        const dso::nanoseconds half =
            dso::nanoseconds(std::lround(dsec.as_underlying_type() / 2));
        dso::datetime<dso::nanoseconds> t = t0;
        t.add_seconds(half);

        // difference in time (aka from tref)
        A(idx, 0) = 1e0;
        const double dt = fit.deltax(t,tref);// t.delta_date(tref).as_mjd();
        A(idx, 1) = dt;
        for (int i = 1; i < poly_order; i++)
          A(idx, i+1) = A(idx, i) * dt;
        // increase matrix index
        ++idx;
      } else if (!(cval % every)) {
        t0 = tend;
      }

      ++cval;
    } // current block within date range
    else if (tl1 >= end) {
      break;
    }
  } // end reading RINEX blocks

  // make sure no error occured while reading the RINEX
  if (error>0) {
    fprintf(stderr,
            "[ERROR] Failed reading \'F\' values off from RINEX %s (traceback: "
            "%s)\n",
            fnrnx, __func__);
    return 1;
  }

  --cval;
  // add any reamaining data in the buffer
  if ((cval % every) != every - 1) {
    const int j = cval % every;
    b(idx) = std::reduce(fs.begin(), fs.end()+j, 0e0) / j;
    
    const auto dsec = tend.delta_sec(t0);
    const dso::nanoseconds half =
        dso::nanoseconds(std::lround(dsec.as_underlying_type() / 2));
    dso::datetime<dso::nanoseconds> t = t0;
    t.add_seconds(half);
    
    // difference in time (aka from tref)
    A(idx,0) = 1e0;
    const double dt = t.delta_date(tref).as_mjd();
    A(idx,1) = dt;
    for (int i = 1; i < poly_order; i++)
      A(idx, i + 1) = A(idx, i) * dt;

    // increase matrix index
    ++idx;
  }

  return 0;
}

int dso::fit_relative_frequency_offset(
    const std::vector<const char *> &fns,
    dso::PolynomialModel<Datetime> &fit,
    bool use_tai,
    int every,
    const Datetime &start,
    const Datetime &end) noexcept 
{
  assert(fit.order() >= 1);

  // Initial guess for the number of RFO values to be read for every RINEX
  // RINEX are (most probably) daily, hence they span 86400 seconds. For every
  // 10sec we have two F values, hence we will have a total of: 8640 * 2 = 
  constexpr const int num_vals_daily = 8640 * 2;

  // we will be taking one F value every N, hence we are going to have a total 
  // of 8640 * 2 / every F data points per day (aka per RINEX)
  const int NDay = std::ceil(num_vals_daily / every);

  // how many days do we want to use data from? If start and end are not set, 
  // then just count the number of RINEX files
  float num_days = fns.size();
  if (start != dso::datetime<dso::nanoseconds>::min() &&
      end != dso::datetime<dso::nanoseconds>::max()) {
    num_days = end.as_mjd() - start.as_mjd();
  }

  // So, this is our educated guess about the number of points we will be 
  // using, for the whole time-span
  long num_pts = NDay * num_days + 10 /* just to be sure ... */;

  // how many columns in design matrix ?
  const int M = fit.order() + 1;

  // allocate matrices (A,b)
  Eigen::VectorXd b(num_pts);
  Eigen::MatrixXd A(num_pts, M);

  // read through all RINEX's in list, and parse data to the b and A matrices.
  // The reference date will be set from the first RINEX
  unsigned index = 0;
  dso::datetime<dso::nanoseconds> tref = dso::datetime<dso::nanoseconds>::min();
  for (const auto sptr : fns) {
    if (parse_rinex_rfos(sptr, use_tai, every, tref, fit, b, A, index,
                         start, end)) {
      fprintf(stderr,
              "[ERROR] Failed fitting Receiver Frequency Offset values! failed "
              "at RINEX file %s (traceback: %s)\n",
              sptr, __func__);
      return 1;
    }
  }

  // OK! all RINEX files parsed. We now need to resize the matrices, given
  // their true sizes
  b.conservativeResize(index, Eigen::NoChange);
  A.conservativeResize(index, Eigen::NoChange);

  printf("Solving LS with b=%ldx%ld and A=%ldx%ld\n", b.rows(), b.cols(), A.rows(), A.cols());

  // Solve the LS problem
  const auto y = (A.transpose() * A).ldlt().solve(A.transpose() * b);

  // assign to the (returned) model
  fit.xref = tref;
  for (int i=0; i<fit.order()+1; i++)
    fit.cf[i] = y(i);

  return 0;
}
