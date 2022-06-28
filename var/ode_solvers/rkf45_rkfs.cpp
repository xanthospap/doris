#include "rkf45.hpp"
#include <cmath>
#include <limits>
#include <cstdio>

int dso::RKF45::rkfs(double &t, double tout, Eigen::VectorXd &y) noexcept {
  
  constexpr const double remin = 1e-12;
  constexpr const int maxnfe = 3000;
  constexpr const double twoeps = 2e0 * std::numeric_limits<double>::epsilon();
  constexpr const double u26 = 13e0 * twoeps;

  // check input parameters for error
  if (neqn < 1 || (iflag < -2 || iflag > 8) || (relerr < 0e0 || abserr < 0e0) ||
      (t == tout && kflag != 3)) {
    iflag = 8;
    return 1;
  }

  int mflag = std::abs(iflag);

  // check for first call, aka iflag == -1 or +1
  if (mflag != 1) {
    // not the first call ....

    // check continuation possibilities
    if (iflag == +2 || iflag == -2) {
      // integration has reached tout
      if ((kflag == 5 && abserr == 0e0) ||
          (kflag == 6 && relerr <= savre && abserr <= savae)) {
        // integration cannot be continued since user did not respond to
        // the instructions pertaining to iflag=5,6,7 or 8
        return 1;
      }
      if (kflag == 3 || !init) {
        iflag = jflag; // 45:
        if (kflag == 3)
          mflag = std::abs(iflag);
      } else if (kflag == 4) {
        // reset function evaluation counter (40:)
        nfe = 0; // 40:
        if (mflag != 2) {
          // reset flag value from previous call (45:)
          iflag = jflag; // 45:
          if (kflag == 3)
            mflag = std::abs(iflag);
        }
      }
    } else {
      // iflag = 3,4,5,6,7 or 8
      if (iflag == 3) {
        // reset flag value from previous call (45:)
        iflag = jflag; // 45:
        if (kflag == 3)
          mflag = std::abs(iflag);
      } else if (iflag == 4) {
        // reset function evaluation counter (40:)
        nfe = 0; // 40:
        if (mflag != 2) {
          // reset flag value from previous call (45:)
          iflag = jflag; // 45:
          if (kflag == 3)
            mflag = std::abs(iflag);
        }
      } else if (iflag == 5 && abserr > 0e0) {
        // reset flag value from previous call (45:)
        iflag = jflag; // 45:
        if (kflag == 3)
          mflag = std::abs(iflag);
      } else {
        // integration cannot be continued since user did not respond to
        // the instructions pertaining to iflag=5,6,7 or 8
        // 30: stop
        return 1;
      }
    }
  }

  // Exiting to 50: ...

  // save input iflag and set continuation flag value for subsequent
  // input checking (50:)
  jflag = iflag;
  kflag = 0;

  // save relerr and abserr for checking input on subsequent calls
  savre = relerr;
  savae = abserr;

  // restrict relative error tolerance to be at least as large as
  // 2*eps+remin to avoid limiting precision difficulties arising
  // from impossible accuracy requests
  const double rer = twoeps + remin;
  if (!(relerr >= rer)) {
    // relative error tolerance too small
    relerr = rer;
    iflag = 3;
    kflag = 3;
    return 0;
  }

  // continue (55:)
  double dt = tout - t;
 
  //  initialization (60:)
  //  --------------------------------------------------------------------
  //  set initialization completion indicator,init
  //  set indicator for too many output points,kop
  //  evaluate initial derivatives
  //  set counter for function evaluations,nfe
  //  evaluate initial derivatives
  //  set counter for function evaluations,nfe
  //  estimate starting stepsize
  init = 0;
  kop = 0;
  double a = t;
  f(a, y, yp, params);
  nfe = 1;
  if (t == tout) {
    iflag = 2;
    return 0;
  }

  // break point 65 (65:)
  init = 1;
  h = std::abs(dt);
  double toln = 0e0;
  for (int i = 0; i < neqn; i++) {
    const double tol = relerr * std::abs(y(i)) + abserr;
    if (tol > 0e0) {
      toln = tol;
      const double ypk = std::abs(yp(i));
      if (ypk * std::pow(h, 5e0) > tol) {
        h = (tol / ypk) * (tol / ypk);
      }
    }
  }
  if (toln <= 0e0) {
    h = 0e0;
  }
  h = std::max(h, u26 * std::max(std::abs(t), std::abs(dt)));
  jflag = std::copysign(2, iflag);

  // set stepsize for integration in the direction from t to tout (80:)
  h = std::copysign(h, dt);

  // test to see if rkf45 is being severely impacted by too many
  // output points
  if (std::abs(h) >= 2e0 * std::abs(dt))
    ++kop;
  if (kop >= 100) {
    // unnecessary frequency of output
    kop = 0;
    iflag = 7;
    return 0;
  }

  // break point 85:
  if (std::abs(dt) <= u26 * std::abs(t)) {
    // if too close to output point,extrapolate and return
    y += dt * yp;
    a = tout;
    f(a, y, yp, params);
    ++nfe;
    // goto 300
    t = tout;
    iflag = 2;
    return 0;
  }

  // initialize output point indicator (95:)
  int output = false;

  // to avoid premature underflow in the error tolerance function,
  // scale the error tolerances
  const double scale = 2e0 / relerr;
  const double ae = scale * abserr;

  // step by step integration (100:)
  do {
    int hfaild = false;

    // set smallest allowable stepsize
    const double hmin = u26 * std::abs(t);

    // adjust stepsize if necessary to hit the output point.
    // look ahead two steps to avoid drastic changes in the stepsize and
    // thus lessen the impact of output points on the code.
    dt = tout - t;
    if (std::abs(dt) >= 2e0 * std::abs(h)) {
      ;
    } else if (std::abs(dt) > std::abs(h)) {
      // goto 150;
      h = 0.5e0 * dt;
    } else {
      // the next successful step will complete the integration to the
      // output point
      output = true;
      h = dt;
    }

    // core integrator for taking a single step (200:)
    // ---------------------------------------------------------------------
    // the tolerances have been scaled to avoid premature underflow in
    // computing the error tolerance function et.
    // to avoid problems with zero crossings,relative error is measured
    // using the average of the magnitudes of the solution at the
    // beginning and end of a step.
    // the error estimate formula has been grouped to control loss of
    // significance.
    // to distinguish the various arguments, h is not permitted
    // to become smaller than 26 units of roundoff in t.
    // practical limits on the change in the stepsize are enforced to
    // smooth the stepsize selection process and to avoid excessive
    // chattering on problems having discontinuities.
    // to prevent unnecessary failures, the code uses 9/10 the stepsize
    // it estimates will succeed.
    // after a step failure, the stepsize is not allowed to increase for
    // the next attempted step. this makes the code more efficient on
    // problems having discontinuities and more effective in general
    // since local extrapolation is being used and extra caution seems
    // warranted.

    // test number of derivative function evaluations.
    // if okay,try to advance the integration from t to t+h (200:)
    double s, esttol;
    do {

      // too much work!
      if (nfe > maxnfe) {
        iflag = 4;
        kflag = 4;
        return 0;
      }

      // advance an approximate solution over one step of length h (220:)
      F.col(0) = fehl(t, h, y, yp);
      nfe += 5;

      // compute and test allowable tolerances versus local error estimates
      // and remove scaling of tolerances. note that relative error is
      // measured with respect to the average of the magnitudes of the
      // solution at the beginning and end of the step.
      double eeoet = 0e0;
      for (int i = 0; i < neqn; i++) {
        const double et = std::abs(y(i)) + std::abs(F(0, i)) + ae;
        if (et <= 0e0) {
          iflag = 5;
          return 0;
        }
        // 240:
        const double ee = std::abs(
            (-2090e0 * yp(i) + (21970e0 * F(i, 2) - 15048e0 * F(i, 3))) +
            (22528e0 * F(i, 1) - 27360e0 * F(i, 4)));
        eeoet = std::max(eeoet, ee / et);
      }

      esttol = std::abs(h) * eeoet * scale / 752400e0;
      if (esttol <= 1e0) {
        // successeful step, break loop at 200: and goto 260;
        break;
      }

      // unsuccessful step
      // reduce the stepsize , try again
      // the decrease is limited to a factor of 1/10
      hfaild = true;
      output = false;
      s = 1e-1;
      if (esttol < 59049e0)
        s = 9e-1 / (esttol * esttol);
      h *= s;

      // if (std::abs(h) > hmin)
      //   goto 200;

      if (std::abs(h) <= hmin) {
        // requested error unattainable at smallest allowable stepsize
        iflag = 6;
        kflag = 6;
        return 0;
      }

    } while (std::abs(h) > hmin);

    // successful step
    // store solution at t+h and evaluate derivatives there (260:)
    t += h;
    // for (int k = 0; k < neqn; k++)
    //   y(k) = F(0,k);
    y = F.col(0);
    a = t;
    f(a, y, yp, params);
    ++nfe;

    // choose next stepsize
    // the increase is limited to a factor of 5
    // if step failure has just occurred, next
    // stepsize is not allowed to increase
    s = 5e0;
    if (esttol > 1.889568e-4)
      s = 9e-1 / (esttol * esttol);
    if (hfaild)
      s = std::min(s, 1e0);
    h = std::copysign(std::max(s * std::abs(h), hmin), h);

    // end of core integrator

    // should we take another step
    if (output) {
      // goto 300;
      t = tout;
      iflag = 2;
      return 0;
    }
    // if (iflag > 0)
    //  goto 100;
  } while (iflag > 0);

  // integration successfully completed

  // one-step mode
  iflag = -2;
  return 0;

  // interval mode (300:)
  // t = tout;
  // iflag = 2;
  // return 0;
}
