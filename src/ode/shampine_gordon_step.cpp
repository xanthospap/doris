// Function prototype for first order differential equations
// void f (double x, const Vector y, Vector yp[])
typedef void (*ODEfun)(
  double        x,     // Independent variable
  const double *y,     // State vector 
  double       *yp,    // Derivative y'=f(x,y)
);

class GSIntegrator {
  // C++ implementation
  // https://github.com/xrf/sg-ode/blob/master/sg_ode/ode.c
  // also
  // https://people.math.sc.edu/Burkardt/cpp_src/ode/ode.html
  private:
    ODEfun f;
    int neqn;
    double g[14],sig[14];
    double phi; // << TODO

    static constexpr const double gstr[/*13*/] = {
        0.5e00,    0.0833e0,  0.0417e0,   0.0264e0,  0.0188e0,
        0.0143e0,  0.0114e0,  0.00936e0,  0.00789e0, 0.00679e0,
        0.00592e0, 0.00524e0, 0.00468e0};
    static constexpr const double two[/*13*/] = {
        2e0,   4e0,   8e0,    16e0,   32e0,   64e0,  128e0,
        256e0, 512e0, 1024e0, 2048e0, 4096e0, 8192e0};

  public:
    int step();
    int interp(double x, const double *y, double xout, double *yout,
                     double *ypout, int neqn, int kold, const double *phi,
                     const double *psi) noexcept;
};

constexpr const double meps = std::numeric_limits<double>::epsilon();

// STEP integrates the system of ODEs one step, from X to X+H.
// 
int dso::GSIntegrator::step(double &x, const double *y, double &h, double &eps,
                            const double *wt, int &start, double &hold,
                            int &crash, int &phase1, int &ns, int &nornd) {
  //                                                                   
  // Begin block 0                                                     
  //                                                                   
  // Check if step size or error tolerance is too small for machine    
  // precision.  If first step, initialize phi array and estimate a    
  // starting step size. If step size is too small, determine an       
  // acceptable one.                                                   
  //                                                                   
  if (std::abs(h) < 4e0 * meps * std::abs(x)) {
    h = std::copysign(4e0 * meps * std::abs(x), h);
    crash = true;
    return;
  }
  crash  = false;

  // If error tolerance is too small, increase it to an
  // acceptable value.
  double round = 0e0;
  for (int i=0; i<neqn; i++) round += (y[i]*y[i]) / (wt[i]*wt[i]);
  round = 2e0 * meps * std::sqrt(round);
  if (0.5e0*eps < round) {
    eps = 2e0 * round* (1e0+4e0*meps);
    crash = true;
    return;
  }
  
  g[1]   = 1e0;
  g[2]   = 0e5;
  sig[1] = 1e0;
  
  if (start) {
    // Initialize. Compute appropriate step size for first step
    f(x,y,yp);
    double sum = 0e0;
    for (int i=0; i<neqn; i++) {
      // TODO
      phi(i,1) = yp[i];
      phi(i,2) = 0e0;
      sum += (yp[i]*yp[i]) / (wt[i]*wt[i]);
    }
    sum = std::sqrt(sum);
    const double absh = std::abs(h);
    if (eps < 16e0 * sum * h *h) absh = 0.25e0 * std::sqrt(eps/sum);
    h = std::copysign(std::max(absh, 4e0*meps*std::abs(x)), h);
    hold = 0e0;
    double hnew = 0e0;
    int k = 1;
    int kold = 0;
    start = false;
    phase1 = true;
    nornd = true;
    if (0.5e0*eps <= 1e2 * round) {
      nornd = false;
      // TODO
      for (int i=0; i<neqn; i++) phi(i,15) = 0e0;
    }
  } // end start
  ifail = 0;

  //
  // Begin block 1
  //
  // compute coefficients of formulas for this step. Avoid computing
  // those quantities not changed when step size is not changed.

  // ns is the number of steps taken with size h, including the current
  // one.  when k < ns, no coefficients change
  if (h != hold) ns = 0;
  if (ns <= kold) ++ns;
  const int nsp1 = ns + 1;
  
  if (k >= ns) {
    // PRE: ns => 1
    // Compute those components of alpha[*],beta[*],psi[*],sig[*]
    // which are changed
    double temp1 = h * (double)ns;
    beta[ns] = 1e0;
    alpha[ns] = 1e0 / (double)ns;
    sig[ns+1] = 1e0;
    if (k>=nsp1) {
      for (int i=ns+1; i<k+1; i++) {
        const double psiim1 = psi[i-1];
        psi[i-1] = temp1;
        beta[i] = beta[i-1] * psi[i-1] / psiim1;
        temp1 = psiim1 + h;
        alpah[i] = h / temp1;
        sig[i+1] = (double)i * alpha[i]*sig[i];
      }
    } // k >= nsp1
    psi[k] = temp1;

    // Compute coefficients g[*]; initialize v[*] and set w[*].
    if (ns>1) {
      // If order was raised, update diagonal part of v[*] 
      if (k>kold) {
        v[k] = 1e0 / (k * (k+1));
        int nsm2 = ns-2;
        for (int j=1; j<ns-1; j++) {
          int i = k-j;
          v[i] -= alpha[j+1] * v[i+1];
        }
      } // k > kold
      // Update V[*] and set W[*] 
      const int limit1 = kp1 - ns;
      const double temp5 = alpha[ns];
      for (int iq=1; iq<limit1+1; iq++) {
        v[iq] -= temp5 * v[iq+1];
        w[iq] = v[iq];
      }
      g[ns+1] = w[1];
    } else { // !(ns>1)
      for (int iq=1; iq<k+1; iq++) {
        const double temp = iq*(iq+1);
        v[iq] = 1e0 / temp;
        w[iq] = v[iq];
      }
    }

    // Compute the g[*] in the work vector w[*]
    if (kp1 >= ns+2) {
      for (int i=ns+2; i<kp1+1; i++) {
        const double limit2 = kp - i;
        const double temp6 = alpha[i-1];
        for (int iq=1; iq<limit2+1; iq++) {
          w[iq] -= temp6 * w[iq+1];
        }
        g[i] = w[1];
      }
    } // kp1 >= ns+2
  } // ks >= ns

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
      // TODO
      for (int l = 0; l < neqn; l++)
        phi(l, i) *= beta[i];
    }
  }
  // Predict solution and differences
  // φ[k+1] <- φ[k] and
  // φ[k] <- 0
  for (int i = 0; i < neqn; i++) {
    // TODO
    phi(i, k+2) = phi(i, k+1);
    phi(i, k+1) = 0e0;
    p[l] = 0e0;
  }
  // 𝐩 ← 𝐩 + g[i] 𝛗[i]
  // 𝛗[i] ← 𝛗[i] + 𝛗[i + 1]
  for (int j = 1; j < k+1; j++) {
    for (int l = 0; l < neqn; l++) {
      p[l] = p[l] + g[i] * phi(l, i);
      phi(l, i) += phi(l, i+1);
    }
  }
  if (nornd) {
    // 𝐩 ← 𝐲 + h 𝐩
    p = y + h * p;
  } else {
    // 𝛕 = h 𝐩 - 𝛗[14]
    // 𝐩 ← 𝐲 + 𝛕
    // 𝛗[15] ← (𝐩 - 𝐲) - 𝛕
    for (int l = 0; l < neqn; l++) {
      const double tau = h * p[l] - phi(l, 15);
      p[l] = y[l] + tau;
      phi(l, 16) = (p[l] - y[l]) - tau;
    }
  }
  xold = x;
  x += h;
  absh = std::abs(h);
  f(x, p, yp);

  // Estimate errors at orders k, k-1, k-2
  // erkm2 = ‖(𝛗[k - 2] + 𝐲′ - 𝛗[0]) / 𝛚‖²
  // erkm1 = ‖(𝛗[k - 1] + 𝐲′ - 𝛗[0]) / 𝛚‖²
  // erk   = ‖(𝐲′ - 𝛗[0]) / 𝛚‖²
  double erkm2 = 0e0;
  double erkm1 = 0e0;
  double erk   = 0e0;
  for (int l = 0; l < neqn; l++) {
    const double temp3 = 1e0 / wt[l];
    const double temp4 = yp[l] - phi(l, 1);
    if (km2 > 0)
      erkm2 += 
              ((phi(l, km1) + temp4) * temp3) * ((phi(l, km1) + temp4) * temp3);
    if (km2 >= 0)
      erkm1 += ((phi(l, k) + temp4) * temp3) * ((phi(l, k) + temp4) * temp3);
    erk += (temp4 * temp3) * (temp4 * temp3);
  }

  if (km2 > 0)
    erkm2 = absh * sig[km1] * gstr[km2] * std::sqrt(erkm2);
  if (km2 >= 0)
    erkm1 = absh * sig[k] * gstr[km1] * std::sqrt(erkm1);

  const double temp5 = absh*std::sqrt(erk);
  err = temp5*(g[k]-g[kp1]);
  erk = temp5*sig[kp1]*gstr[k];
  knew = k;

  // Test if order should be lowered
  if (km2 > 0)
    if (std::max(erkm1, erkm2) <= erk)
      knew = km1;
  if (km2 == 0)
    if (erkm1 <= 0.5e0 * erk)
      knew = km1;

  //
  // If step is successful continue with block 4, otherwise repeat
  // blocks 1 and 2 after executing block 3
  //
  int success = (err<=eps);
  if (!success) {
      //
      // Begin block 3
      //

      // The step is unsuccessful. Restore x, phi[*,*], psi[*]. If
      // 3rd consecutive failure, set order to 1. If step fails more
      // than 3 times, consider an optimal step size. Double error
      // tolerance and return if estimated step size is too small
      // for machine precision.
      //

  }
}
