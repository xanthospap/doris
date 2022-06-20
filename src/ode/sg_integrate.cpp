#include "ode.hpp"

constexpr const double FourEps = 4e0 * std::numeric_limits<double>::epsilon();

dso::EnumBitset<dso::ODEStatusCode>
dso::GSOdeSolver::integrate(double &t, double tout, double *y,
                            int init) noexcept {

  dso::EnumBitset<dso::ODEStatusCode> ebs;

  // Test for improper parameters
  double eps = std::max(relerr, abserr);

  if ((relerr < 0e0) || // Negative relative error bound
      (abserr < 0e0) || // Negative absolute error bound
      (eps <= 0e0) ||   // Both error bounds are non-positive
      ((!init) && (t != told))) {
    return ebs.set(dso::ODEStatusCode::INVPARAM); // Set error code
  }

  // On each call set interval of integration and counter for
  // number of steps. Adjust input error tolerances to define
  // weight vector for subroutine STEP.
  const double del = tout - t;
  const double absdel = std::abs(del);

  const double tend = t + 1e2 * del;
  if (!permitTout)
    tend = tout;

  int nostep = 0;
  int kle4 = 0;
  int stiff = false;
  double releps = relerr / eps;
  double abseps = abserr / eps;

  if ((init) || (!oldPermit) || (delsgn * del <= 0e0)) {
    // On start and restart also set the work variables x and yy(*),
    // store the direction of integration and initialize the step size
    start = true;
    x = t;
    for (int i = 0; i < neqn; i++)
      yy(i) = y[i];
    delsgn = std::copysign(1e0, del);
    h = std::copysign(std::max(FourEps * std::abs(x), std::abs(tout - x)),
                      tout - x);
  }

  while (true) {
    // If already past output point, interpolate solution and return
    if (std::abs(x - t) >= absdel) {
      interpolate(tout, y, ypout);
      t = tout; // Set independent variable
      told = t; // Store independent variable
      oldPermit = permitTout;
      // Normal exit
      return ebs.set(dso::ODEStatusCode::DONE);
    }

    // If cannot go past output point and sufficiently close,
    // extrapolate and return
    if (!permitTout && (std::abs(tout - x) < FourEps * std::abs(x))) {
      h = tout - x;
      f(x, yy, yp);    // Compute derivative yp(x)
      y = yy + h * yp; // Extrapolate vector from x to tout
      t = tout;        // Set independent variable
      told = t;        // Store independent variable
      oldPermit = permitTout;
      // Normal exit
      return ebs.set(dso::ODEStatusCode::DONE);
    }

    // Test for too much work
    if (nostep >= maxnum) {
      // State = DE_NUMSTEPS;          // Too many steps
      if (stiff)
        ebs.set(dso::ODEStatusCode::DE_STIFF); // Stiffness suspected
      y = yy;                                  // Copy last step
      t = x;
      told = t;
      OldPermit = true;
      return ebs.set(dso::ODEStatusCode::NUMSTEPS); // Weak failure exit
    }

    // Limit step size, set weight vector and take a step
    h = std::copysign(std::min(std::abs(h), std::abs(tend - x)), h);
    for (int l = 0; l < neqn; l++)
      wt(l) = releps * std::abs(yy(l)) + abseps;
    step(x, yy, eps, crash);

    // Test for too small tolerances
    if (crash) {
      relerr = eps * releps; // Modify relative and absolute
      abserr = eps * abseps; // accuracy requirements
      y = yy;                // Copy last step
      t = x;
      told = t;
      oldPermit = true;
      return ebs.set(dso::ODEStatusCode::BADACC); // Weak failure exit
    }

    ++nostep; // Count total number of steps

    // Count number of consecutive steps taken with the order of
    // the method being less or equal to four and test for stiffness
    ++kle4;
    if (kold > 4)
      kle4 = 0;
    if (kle4 >= 50)
      stiff = true;

  } // end step loop

  return ebs;
}
