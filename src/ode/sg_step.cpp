#include "ode.hpp"
#ifdef DEBUG
#include <cassert>
#endif

constexpr const double meps = std::numeric_limits<double>::epsilon();
constexpr const double FourEps = 4e0 * std::numeric_limits<double>::epsilon();

// STEP integrates the system of ODEs one step, from X to X+H.
//
int dso::GSOdeSolver::step(double &x, double *y, double &eps,
                           int &crash) noexcept {
  const double p5eps = 0.5e0 * eps;

  //
  // Begin block 0
  //
  // Check if step size or error tolerance is too small for machine
  // precision.  If first step, initialize phi array and estimate a
  // starting step size. If step size is too small, determine an
  // acceptable one.
  //
  if (std::abs(h) < FourEps * std::abs(x)) {
    h = std::copysign(FourEps * std::abs(x), h);
    crash = true;
    return;
  }
  crash = false;

  // If error tolerance is too small, increase it to an
  // acceptable value.
  // round = 2 ε ‖𝐲 / 𝛚‖
  double round = 0e0;
  for (int i = 0; i < neqn; i++)
    round += (y[i] * y[i]) / (wt[i] * wt[i]);
  round = 2e0 * meps * std::sqrt(round);
  if (p5eps < round) {
    eps = 2e0 * round * (1e0 + FourEps);
    crash = true;
    return;
  }

  g[1] = 1e0;
  g[2] = 0e5;
  sig[1] = 1e0;

  if (start) {
    // Initialize. Compute appropriate step size for first step
    f(x, y, yp);
    double sum = 0e0;
    for (int i = 0; i < neqn; i++) {
      phi(i, 1) = yp[i];
      phi(i, 2) = 0e0;
      sum += (yp[i] * yp[i]) / (wt[i] * wt[i]);
    }
    sum = std::sqrt(sum);
    double absh = std::abs(h);
    if (eps < 16e0 * sum * h * h)
      absh = 0.25e0 * std::sqrt(eps / sum);
    h = std::copysign(std::max(absh, 4e0 * meps * std::abs(x)), h);
    hold = 0e0; /* instance member */
    double hnew = 0e0;
    k = 1;         /* instance member */
    kold = 0;      /* instance member */
    start = false; /* instance member */
    phase1 = true; /* instance member */
    nornd = true;  /* instance member */
    if (p5eps <= 1e2 * round) {
      nornd = false;
      // 𝛗[14] ← 𝟎
      for (int i = 0; i < neqn; i++)
        phi(i, 15) = 0e0;
    }
  } // end start
  ifail = 0;

#ifdef DEBUG
  int kp1;
#endif
  //
  // Repeat blocks 1, 2 (and 3) until step is successful
  //
  do {
    kp1 = k + 1;

    //
    // Begin block 1
    //
    // compute coefficients of formulas for this step. Avoid computing
    // those quantities not changed when step size is not changed.

    // ns is the number of steps taken with size h, including the current
    // one.  when k < ns, no coefficients change
    if (h != hold)
      ns = 0;
    if (ns <= kold)
      ++ns;
    const int nsp1 = ns + 1;

    if (k >= ns) {
      // PRE: ns => 1
      // Compute those components of alpha[*],beta[*],psi[*],sig[*]
      // which are changed
      double temp1 = h * (double)ns;
      beta(ns) = 1e0;
      alpha(ns) = 1e0 / (double)ns;
      sig(ns + 1) = 1e0;
      if (k >= ns + 1) {
        for (int i = ns + 1; i < k + 1; i++) {
          const double temp2 = psi(i - 1);
          psi(i - 1) = temp1;
          beta(i) = beta(i - 1) * psi(i - 1) / temp2;
          temp1 = temp2 + h;
          alpah(i) = h / temp1;
          sig(i + 1) = (double)i * alpha(i) * sig(i);
        }
      } // k >= nsp1
      psi(k) = temp1;

      // Compute coefficients g[*]; initialize v[*] and set w[*].
      if (ns > 1) {
        // If order was raised, update diagonal part of v[*]
        if (k > kold) {
          v(k) = 1e0 / (k * (k + 1));
          for (int j = 1; j < ns - 1; j++) {
            int i = k - j;
            v(i) -= alpha(j + 1) * v(i + 1);
          }
        } // k > kold
        // Update V[*] and set W[*]
        const int limit1 = k + 1 - ns;
        for (int iq = 1; iq < limit1 + 1; iq++) {
          v(iq) -= alpha(ns) * v(iq + 1);
          w(iq) = v(iq);
        }
        g(ns + 1) = w(1);
      } else { // !(ns>1)
        for (int iq = 1; iq < k + 1; iq++) {
          const double temp = iq * (iq + 1);
          v(iq) = 1e0 / temp;
          w(iq) = v(iq);
        }
      }

      // Compute the g[*] in the work vector w[*]
      if (k + 1 >= ns + 2) {
        for (int i = ns + 2; i < k + 2; i++) {
          const int limit2 = k + 2 - i;
          for (int iq = 1; iq < limit2 + 1; iq++) {
            w(iq) -= alpha(i - 1) * w(iq + 1);
          }
          g(i) = w(1);
        }
      } // kp1 >= ns+2
    }   // ks >= ns

    //
    // Begin block 2
    //
    // Predict a solution p[*], evaluate derivatives using predicted
    // solution, estimate local error at order k and errors at orders
    // k, k-1, k-2 as if constant step size were used.
    //

    // Change phi to phi star
    if (k >= ns + 1) {
      // φ[i] <- φ[i] β[i]
      for (int i = ns + 1; i < k + 1; i++) {
        for (int l = 0; l < neqn; l++) {
          phi(l, i) *= beta[i];
        }
      }
    }

    // Predict solution and differences
    // φ[k+1] <- φ[k] and
    // φ[k] <- 0
    for (int i = 0; i < neqn; i++) {
      phi(i, k + 2) = phi(i, k + 1);
      phi(i, k + 1) = 0e0;
    }

    // 𝐩 ← 𝐩 + g[i] 𝛗[i]
    // 𝛗[i] ← 𝛗[i] + 𝛗[i + 1]
    for (int i = 0; i < neqn; i++)
      p(i) = 0e0;
    for (int j = 1; j < k + 1; j++) {
      const int i = k + 1 - j;
      for (int l = 0; l < neqn; l++) {
        p(l) += g(i) * phi(l, i);
        phi(l, i) += phi(l, i + 1);
      }
    }
    if (nornd) {
      // 𝐩 ← 𝐲 + h 𝐩
      for (int i = 0; i < ; i++)
        p(i) = y[i] + h * p(i);
    } else {
      // 𝛕 = h 𝐩 - 𝛗[14]
      // 𝐩 ← 𝐲 + 𝛕
      // 𝛗[15] ← (𝐩 - 𝐲) - 𝛕
      for (int l = 0; l < neqn; l++) {
        const double tau = h * p(l) - phi(l, 15);
        p(l) = y[l] + tau;
        phi(l, 16) = (p(l) - y[l]) - tau;
      }
    }
    double xold = x;
    x += h;
    absh = std::abs(h);
    f(x, p, yp);

    // Estimate errors at orders k, k-1, k-2
    // erkm2 = ‖(𝛗[k - 2] + 𝐲′ - 𝛗[0]) / 𝛚‖²
    // erkm1 = ‖(𝛗[k - 1] + 𝐲′ - 𝛗[0]) / 𝛚‖²
    // erk   = ‖(𝐲′ - 𝛗[0]) / 𝛚‖²
    double erkm2 = 0e0;
    double erkm1 = 0e0;
    double erk = 0e0;
    for (int l = 0; l < neqn; l++) {
      const double temp3 = 1e0 / wt(l);
      const double temp4 = yp(l) - phi(l, 1);
      if (k > 2)
        erkm2 += ((phi(l, k - 1) + temp4) * temp3) *
                 ((phi(l, k - 1) + temp4) * temp3);
      if (k >= 2)
        erkm1 += ((phi(l, k) + temp4) * temp3) * ((phi(l, k) + temp4) * temp3);
      erk += (temp4 * temp3) * (temp4 * temp3);
    }

    if (k > 2)
      erkm2 = absh * sig(k - 1) * gstr(k - 2) * std::sqrt(erkm2);
    if (k >= 2)
      erkm1 = absh * sig(k) * gstr(k - 1) * std::sqrt(erkm1);

    const double temp5 = absh * std::sqrt(erk);
    const double err = temp5 * (g(k) - g(k + 1));
    const double erk = temp5 * sig(k + 1) * gstr(k);
    int knew = k;

    // Test if order should be lowered
    if (k > 2)
      if (std::max(erkm1, erkm2) <= erk)
        knew = km1;
    if (k == 2)
      if (erkm1 <= 0.5e0 * erk)
        knew = km1;

    //
    // If step is successful continue with block 4, otherwise repeat
    // blocks 1 and 2 after executing block 3
    //
    int success = (err <= eps);
    if (!success) {
      //
      // Begin block 3
      //
      // Note that this step may change the value of k

      // The step is unsuccessful. Restore x, phi[*,*], psi[*]. If
      // 3rd consecutive failure, set order to 1. If step fails more
      // than 3 times, consider an optimal step size. Double error
      // tolerance and return if estimated step size is too small
      // for machine precision.

      // Restore x, phi[*,*] and psi[*]
      phase1 = false;
      x = xold;
      for (int i = 1; i < k + 1; i++) {
        for (int l = 0; l < neqn; l++) {
          phi(l, i) = (phi(l, i) - phi(l, i + 1)) / beta(i);
        }
      }

      if (k >= 2)
        for (int i = 2; i < k + 1; i++)
          psi(i - 1) = psi(i) - h;

      // On third failure, set order to one.
      // Thereafter, use optimal step size
      ++ifail;
      double temp2 = 0.5e0;
      if (ifail >= 3) {
        knew = 1;
        if (ifail > 3 && p5eps < 0.25e0 * erk) {
          temp2 = std::sqrt(p5eps / erk);
        }
      }

      // WARNING! set new k
      h = temp2 * h;
      k = knew;
      if (std::abs(h) < FourEps * std::abs(x)) {
        crash = true;
        h = std::copysign(FourEps * std::abs(x), h);
        eps += eps;
        return;
      }

      //
      // End block 3, return to start of block 1
      //
    } // if (!success)

  } while (!success);

  //
  // Begin block 4
  //
  // The step is successful. Correct the predicted solution, evaluate
  // the derivatives using the corrected solution and update the
  // differences. Determine best order and step size for next step.
  //
  kold = k;
  hold = h;

#ifdef DEBUG
  assert(kp1 == k + 1);
#endif

  // Correct and evaluate
  const double hgk = h * g(k + 1);
  if (nornd) {
    // 𝐲 ← 𝐩 + h g[k] (𝐲′ - 𝛗[0])
    for (int l = 0; l < neqn; l++)
      y[l] = p(l) + hgk * (yp(l) - phi(l, 1));
  } else {
    // 𝛒 = h g[k] (𝐲′ - 𝛗[0]) - 𝛗[15]
    // 𝐲 ← 𝐩 + 𝛒
    // 𝛗[14] ← (𝐲 - 𝐩) - 𝛒
    for (int l = 0; l < neqn; l++) {
      const double rho = hgk * (yp(l) - phi(l, 1)) - phi(l, 16);
      y[l] = p(l) + rho;
      phi(l, 15) = (y[l] - p(l)) - rho;
    }
  }
  f(x, y, yp);

  // update differences for next step
  // 𝛗[k] ← 𝐲′ - 𝛗[0]
  // 𝛗[k + 1] ← 𝛗[k] - 𝛗[k + 1]
  for (int l = 0; l < neqn; l++) {
    phi(l, k + 1) = yp(l) - phi(l, 1);
    phi(l, k + 2) = phi(l, k + 1) - phi(l, k + 2);
  }
  // ∀ i ∈ [0, k).  𝛗[i] ← 𝛗[i] + 𝛗[k]
  for (int i = 1; i < k + 1; i++) {
    for (int l = 0; l < neqn; l++) {
      phi(l, i) += phi(l, k + 1);
    }
  }

  // Estimate error at order k+1 unless
  // - in first phase when always raise order,
  // - already decided to lower order,
  // - step size not constant so estimate unreliable
  erkp1 = 0e0;
  if ((knew == k - 1) || (k == kmax))
    phase1 = false;

  if (phase1) {
    ++k;
    erk = erkp1;
  } else {
    if (knew == k - 1) {
      // lower order
      --k;
      erk = erkm1;
    } else {
      if (k + 1 <= ns) {
        // erkp1 = ‖𝛗[k + 1] / 𝛚‖²
        for (int l = 0; l < neqn; l++) {
          erkp1 += (phi(l, k + 2) / wt(l)) * (phi(l, k + 2) / wt(l));
        }
        erkp1 = absh * gstr(k + 1) * std::sqrt(erkp1);

        // Using estimated error at order k+1, determine
        // appropriate order for next step
        if (k > 1) {
          if (erkm1 <= std::min(erk, erkp1)) {
            // lower order
            --k;
            erk = erkm1;
          } else {
            if (erkp1 < erk && k != kmax) {
              // raise order
              ++k;
              erk = erkp1;
            }
          }
        } else {
          if (erkp1 < 0e5 * erk) {
            // raise order
            // Here erkp1 < erk < max(erkm1,ermk2) else
            // order would have been lowered in block 2.
            // Thus order is to be raised
            ++k;
            erk = erkp1;
          }
        }
      } // if (kp1<=ns)
    }   // if (knew==km1)
  }     // if (phase1)

  // With new order determine appropriate step size for next step
  if (phase1 || (p5eps >= erk * two[k + 1])) {
    hnew = 2e0 * h;
  } else {
    if (p5eps < erk) {
      const double r = std::pow(p5eps / erk, 1e0 / (k + 1));
      hnew = absh * std::max(0.5e0, std::min(0.9e0, r));
      hnew = std::copysign(std::max(hnew, FourEps * std::abs(x)), h);
    } else {
      hnew = h;
    }
  }

  h = hnew;

  return;
}
