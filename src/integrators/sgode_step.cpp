#include "sgode.hpp"
#include <limits>
#ifdef DEBUG
#include <cassert>
#endif

/* Accoring to Shampine & Gordon, step is the basic Adams integrator
 * Subroutine STEP solves a system of first order ordinary differential
 * equations:
 *                  y'(x) = f(x, y(x))
 *                       y(a) = y0
 * using Adams methods.
 * We anticipate most users will use this routine indirectly through the 
 * driver DE.
 * The code STEP advances the solution of the differential equation one step 
 * and returns control to the calling program. Because of this and because of 
 * the information returned, it allows the integration to be monitored closely. 
 * It also allows great flexibility in the way the local error is controlled.
 * As with DE, the output of one call to STEP is nearly all the input for the 
 * next call so that STEP is only a little more trouble to use than DE.
 *
 * Since STEP advances the solution one step at a time, it must be called 
 * repeatedly to integrate to a specified output point. And so, the typical 
 * situation is that a successful step has just been completed and the user 
 * wishes to take another step. On return from a successful step, X represents 
 * how far the integration has progressed, Y is the vector of solution components
at X, YP is the vector of first derivatives of the solution at X, H is
estimated by the code to be the largest step the code can take and still pass
the error test, START = .FALSE., and CRASH = .FALSE.. The parameter
EPS is the local error tolerance and the vector WT specifies the error
criterion; we shall go into them later. The user may change EPS and WT
as he wishes. For some kinds of error control it is necessary to change WT
at each step. The user may alter H but ordinarily he should not. The code
does a very good job of selecting H and the user should not override the
choice without a good reason. It is most efficient to step past any desired
output point and to use INTRP to obtain the solution and its derivative
there. Sometimes it is not possible to integrate beyond a given point, as the
equation may have a discontinuity or be undefined there, and it is necessary
for the user to alter H so as to land at the desired point. The code will
attempt to step to X + H when it is called. If this is not possible, it will
reduce the step size and continue trying to complete a step until it either
succeeds or realizes that the task is impossible. Thus while the user does
not know how far the code will advance, he does know that it will not go
beyond X + H. By manipulating H the user can arrange to stop at a given
point. We emphasize that it is much more efficient and more accurate to
step past output points and then use INTRP when this is possible. In
summary, to continue after a successful step, the only quantity the user
commonly changes before calling the code again is WT.
 */

constexpr const double umach = std::numeric_limits<double>::epsilon();
constexpr const double twou = 2e0 * umach;
constexpr const double fouru = 4e0 * umach;
/* constexpr const double two[] = {2e0,    4e0,    8e0,   16e0,  32e0,
                                64e0,   128e0,  256e0, 512e0, 1024e0,
                                2048e0, 4096e0, 8192e0}; */
constexpr const double gstr[] = {
    5e-1,      0.0833e0,  0.0417e0,  0.0264e0,  0.0188e0,  0.0143e0, 0.0114e0,
    0.00936e0, 0.00789e0, 0.00679e0, 0.00592e0, 0.00524e0, 0.00468e0};

// input vector y0 is yy
// eps -- local error tolerance
int dso::SGOde::step(double &eps) noexcept {

  // ***     begin block 0     ***
  // check if step size or error tolerance is too small for machine
  // precision.  if first step, initialize phi array and estimate a
  // starting step size.

  // if step size is too small, determine an acceptable one
  if (std::abs(h) < fouru * std::abs(tc)) {
    // set step size for next call/step
    h = std::copysign(fouru * std::abs(tc), h);
    // return signaling a 'crash'
    return 1;
  }

  const double p5eps = 5e-1 * eps;

  // if error tolerance is too small, increase it to an acceptable value
  const double round = twou * (yy().array() / wt().array()).matrix().norm();
  if (p5eps < round) {
    // increase error tolerance
    eps = 2e0 * round * (1e0 + fouru);
    // return signaling a 'crash'
    return 2;
  }

  g(0) = 1e0;
  g(1) = 5e-1;
  sig(0) = 1e0;

  double absh;
  if (iflag == IFLAG::RESTART) {
    // initialize.  compute appropriate step size for first step
    f(tc, yy(), yp(), *params);
    // φ[0] = y'
    Phi.col(0) = yp();
    // φ[1] = 0
    Phi.col(1).setZero();
    // sum = | y'/w |
    const double sum = (yp().array() / wt().array()).matrix().norm();
    absh = std::abs(h);
    if (eps < 16e0 * sum * h * h)
      absh = 0.25e0 * std::sqrt(eps / sum);
    h = std::copysign(std::max(absh, fouru * std::abs(tc)), h);
    hold = 0e0;
    k = 1;
    kold = 0;
    iflag = IFLAG::UNDEFINED;
    phase1 = true; // TODO WTF is this ?
    nornd = true;  // TODO wtf is this ?
    if (p5eps <= 100e0 * round) {
      nornd = false;
      // φ[14] = 0
      Phi.col(14).setZero();
    }
  }
  // Number of failed attempts to procced one step; see Block 3
  int ifail = 0;
  // ***     end block 0     ***

  //
  // Repeat blocks 1, 2 (and 3) until step is successful
  //
  int knew;
  double erkm2, erkm1, erk, err;
  while (true) {
    //  ***     begin block 1     ***

    // ns is the number of steps taken with size h, including the current
    // one. When k<ns, no coefficients change
    if (h != hold) {
      // new step size, reset ns
      ns = 0;
    }
    if (ns <= kold) {
      // number of steps less that current order, increase step counter
      ++ns;
    }

    // k >= ns, means we have to compute coefficients. If we do not meet this
    // condition (aka k >= ns) we exit the block!
    // compute those components of alpha(*),beta(*),psi(*),sig(*) which
    // are changed
    if (k >= ns) {
      // β_{ns}(n+1) = 1
      beta(ns - 1) = 1e0;
      // α_{ns}(n+1) = 1 / ns
      alpha(ns - 1) = 1e0 / ns;
      double tmp1 = h * ns;
      sig(ns) = 1e0;
      for (int i = ns; i < k; i++) {
        // ψ_{i-1}(n)
        const double t2 = psi(i - 1);
        // ψ_{i}(n+1) = i*h -- overwrite value of ψ_{i}(n)
        psi(i - 1) = tmp1;
        // β_i (n+1) = β_{i-1}(n+1) ψ_{i-1}(n+1) / ψ_{i-1}(n)
        beta(i) = beta(i - 1) * psi(i - 1) / t2;
        tmp1 = t2 + h;
        // α_{i}(n+1) = h / ψ_{i}(n+1)
        alpha(i) = h / tmp1;
        sig(i + 1) = i * alpha(i) * sig(i);
      }
      psi(k - 1) = tmp1;

      // compute coefficients g(*)

      // initialize v(*) and set w(*). g(2) is set in data statement
      if (ns <= 1) {
        for (int iq = 0; iq < k; iq++) {
          v(iq) = 1e0 / ((iq + 1e0) * (iq + 2e0));
        }
        std::memcpy(w(), v(), sizeof(double) * k);
      } else { // ns > 1
        // if order was raised, update diagonal part of v(*)
        if (k > kold) {
          v(k - 1) = 1e0 / (k * (k+1));
          for (int j = 1; j < ns - 1; j++) {
            const int i = k - j;
            v(i - 1) -= alpha(j) * v(i);
          }
        }

        // update v(*) and set w(*)
        {
          const double tmp5 = alpha(ns - 1);
          for (int iq = 0; iq < k - ns; iq++) {
            v(iq) -= tmp5 * v(iq + 1);
          }
          std::memcpy(w(), v(), sizeof(double) * (k - ns));
          g(ns) = w(0);
        }
      }

      // compute the g(*) in the work vector w(*)
      for (int i = ns; i < k; i++) {
        for (int iq = 0; iq < k - 1; iq++) {
          w(iq) -= alpha(i) * w(iq + 1);
        }
        g(i + 1) = w(0);
      }
    } // if (k >= ns)
    // ***     end block 1     **

    // ***     begin block 2     ***
    // predict a solution p(*), evaluate derivatives using predicted
    // solution, estimate local error at order k and errors at orders k,
    // k-1, k-2 as if constant step size were used.

    // Change phi to phi star φ[i] = β[i] * φ[i]
    for (int i = ns; i < k; i++) {
      const double b = beta(i);
      Phi.col(i) *= b;
    }

    // Predict solution and differences
    Phi.col(k + 1) = Phi.col(k); // φ[k+1] = φ[k]
    Phi.col(k).setZero();        // φ[k] = 0
    p().setZero();               // p = φ * g

    // p = p + g[i] * φ[i]
    for (int i = k - 1; i >= 0; i--)
      p() += g(i) * Phi.col(i);
    // φ[i] = φ[i] + φ[i+1]
    for (int i = k - 1; i >= 0; i--)
      Phi.col(i) += Phi.col(i + 1);

    if (!nornd) {
      const auto tau = h * p() - Phi.col(14); // τ = h * p - φ[14]
      p() = yy() + tau;                       // p = y + τ
      Phi.col(15) = (p() - yy()) - tau;       // φ[15] = (p-y) - τ
    } else {
      // p = y + h * p
      p() = yy() + h * p();
    }

    {
      const double tcold = tc;
      tc += h;
      absh = std::abs(h);
      f(tc, p(), yp(), *params);

      // Estimate errors at orders k,k-1,k-2
      erkm2 = erkm1 = erk = 0e0;
      {
        // array y' - φ[0] , temporary
        const Eigen::ArrayXd ypmf0 = yp().array() - Phi.col(0).array();
        // erk   = | (y' - φ[0]) / wt |^2
        Eigen::ArrayXd tmpar = ypmf0 / wt().array();
        erk = tmpar.matrix().squaredNorm();
        // erkm2 = | (φ[k-2] + y' - φ[0]) / wt |^2
        if (k > 2) {
          tmpar = (Phi.col(k - 2).array() + ypmf0) / wt().array();
          erkm2 = absh * sig(k - 2) * gstr[k - 3] * tmpar.matrix().norm();
        }
        // erkm1 = | (φ[k-1] + y' - φ[0]) / wt |^2
        if (k >= 2) {
          tmpar = (Phi.col(k - 1).array() + ypmf0) / wt().array();
          erkm1 = absh * sig(k - 1) * gstr[k - 2] * tmpar.matrix().norm();
        }
      }
      {
        const double t5 = absh * std::sqrt(erk);
        err = t5 * (g(k - 1) - g(k));
        erk = t5 * sig(k) * gstr[k - 1];
        knew = k;
      }

      // test if order should be lowered
      if (k == 2) {
        if (erkm1 <= erk * 0.5e0) {
          knew = k - 1;
        }
      } else if (k > 2) {
        if (std::max(erkm1, erkm2) <= erk) {
          knew = k - 1;
        }
      }

      if (err <= eps) {
        fprintf(stderr, "[DEBUG] Step success! order=%d step=%.12f\n", k, h);
        break;
      }
      // ***     end block 2     ***

      //
      // ***     begin block 3     ***
      // the step is unsuccessful.  restore  x, phi(*,*), psi(*) .
      // if third consecutive failure, set order to one.  if step fails more
      // than three times, consider an optimal step size.  double error
      // tolerance and return if estimated step size is too small for machine
      // precision.
      // ***

      // restore x, φ and ψ
      phase1 = false;
      tc = tcold;
    }
    for (int i = 0; i < k; i++) {
      const double t1 = 1e0 / beta(i);
      Phi.col(i) = t1 * (Phi.col(i) - Phi.col(i + 1));
    }
    for (int i = 1; i < k; i++) {
      psi(i - 1) = psi(i) - h;
    }

    // On third failure, set order to one. Thereafter, use optimal step
    // size
    ++ifail;
    {
      double tmp2 = .5e0;
      if (ifail >= 3) {
        fprintf(stderr, "[DEBUG] Third failed attempt to step forward!\n");
        if (ifail != 3 && p5eps < erk * .25e0) {
          tmp2 = std::sqrt(p5eps / erk);
        }
        knew = 1;
      }
      h = tmp2 * h;
    }
    k = knew;

    if (std::abs(h) < fouru * std::abs(tc)) {
      h = std::copysign(fouru * std::abs(tc), h);
      eps += eps;
      return 3;
    }
    // ***     end block 3     ***
  }

  //
  // ***    begin block 4     ***
  //
  // the step is successful. correct the predicted solution, evaluate
  // the derivatives using the corrected solution and update the
  // differences.  determine best order and step size for next step.
  //
  kold = k;
  hold = h;

  // correct and evaluate
  {
    const double fac = h * g(k);
    if (!nornd) {
      // ρ = h g[k] (y' -φ[0]) - φ[15]
      const auto rho =
          fac * (/*ArraysNeqn.col(3)*/ yp() - Phi.col(0)) - Phi.col(15);
      // y = p + ρ
      yy() = p() + rho;
      // φ[14] = (y-p) - ρ
      Phi.col(14) = (yy() - p()) - rho;
    } else {
      // y = p + h * g[k] (y' - φ[0])
      yy() = p() + fac * (yp() - Phi.col(0));
    }
  }
  f(tc, yy(), yp(), *params);

  // update differences for next step
  // φ[k] = y' - φ[0]
  Phi.col(k) = yp() - Phi.col(0);
  // φ[k+1] = φ[k] - φ[k+1]
  Phi.col(k + 1) = Phi.col(k) - Phi.col(k + 1);
  for (int i = 0; i < k; i++)
    Phi.col(i) += Phi.col(k);

  // Estimate error at order k+1 unless:
  //  * in first phase when always raise order,
  //  * already decided to lower order,
  //  * step size not constant so estimate unreliable
  double erkp1 = 0e0;
  if (knew == k - 1 || k == 12)
    phase1 = false;

  if (phase1) {
    ++k;
    erk = erkp1;
  } else if (knew == k - 1) {
    --k;
    erk = erkm1;
  } else if (k < ns) {
    // erkp1 = |φ[k+1] / wt|^2
    erkp1 = (Phi.col(k + 1).array() / wt().array()).matrix().squaredNorm();
    // using estimated error at order k+1, determine appropriate order
    // for next step
    if (k > 1) {
      if (erkm1 <= std::min(erk, erkp1)) {
        --k;
        erk = erkm1;
      } else if (erkp1 < erk && k != 12) {
        ++k;
        erk = erkp1;
      }
    } else if (erkp1 < erk * .5e0) {
      ++k;
      erk = erkp1;
    }
  }

  // with new order determine appropriate step size for next step
  {
    double hnew = 2e0 * h;
    if (!phase1 && p5eps < erk * std::pow(2, k + 1)) {
      hnew = h;
      if (p5eps >= erk)
        return 0;
      const double r = std::pow(p5eps / erk, 1e0 / (k + 1));
      hnew = absh * std::max(.5e0, std::min(.9e0, r));
      hnew = std::copysign(std::max(hnew, fouru * std::abs(tc)), h);
    }
    h = hnew;
  }

  return 0;
}
