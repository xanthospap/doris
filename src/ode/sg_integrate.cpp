#include "ode.hpp"
#include <cmath>
#include <limits>

constexpr const double FourEps = 4e0 * std::numeric_limits<double>::epsilon();

dso::EnumBitset<dso::ODEStatusCode>
dso::GSOdeSolver::integrate(double &t, double tout, double *y) noexcept {

  dso::EnumBitset<dso::ODEStatusCode> ebs;

  // Test for improper parameters
  double eps = std::max(this->relerr, this->abserr);

  if ((this->relerr < 0e0) || // Negative relative error bound
      (this->abserr < 0e0) || // Negative absolute error bound
      (eps <= 0e0) ||         // Both error bounds are non-positive
      ((!this->init) && (t != told))) {
    return ebs.set(dso::ODEStatusCode::INVPARAM); // Set error code
  }

  // On each call set interval of integration and counter for
  // number of steps. Adjust input error tolerances to define
  // weight vector for subroutine STEP.
  const double del = tout - t;
  const double absdel = std::abs(del);

  const double tend = (permitTout) ? (t + 1e2 * del) : (tout);

  int nostep = 0;
  int kle4 = 0;
  int stiff = false;
  double releps = this->relerr / eps;
  double abseps = this->abserr / eps;

  if ((this->init) || (!this->oldPermit) || (this->delsgn * del <= 0e0)) {
    // On start and restart also set the work variables x and yy(*),
    // store the direction of integration and initialize the step size
    this->start = true;
    this->x = t;
    for (int i = 0; i < neqn; i++)
      this->yy(i) = y[i];
    this->delsgn = std::copysign(1e0, del);
    this->h = std::copysign(
        std::max(FourEps * std::abs(this->x), std::abs(tout - this->x)),
        tout - this->x);
  }

  while (true) {
    // If already past output point, interpolate solution and return
    if (std::abs(this->x - t) >= absdel) {
      interpolate(tout, y, this->ypout());
      t = tout; // Set independent variable
      told = t; // Store independent variable
      this->oldPermit = this->permitTout;
      // Normal exit
      return ebs.set(dso::ODEStatusCode::DONE);
    }

    // If cannot go past output point and sufficiently close,
    // extrapolate and return
    if (!this->permitTout &&
        (std::abs(tout - this->x) < FourEps * std::abs(this->x))) {
      h = tout - this->x;
      f(this->x, this->yy(), this->yp()); // Compute derivative yp(x)
      for (int i = 0; i < neqn; i++)
        y[i] = yy(i) + h * yp(i); // Extrapolate vector from x to tout
      t = tout;                   // Set independent variable
      told = t;                   // Store independent variable
      this->oldPermit = this->permitTout;
      // Normal exit
      return ebs.set(dso::ODEStatusCode::DONE);
    }

    // Test for too much work
    if (nostep >= this->maxnum) {
      // Too many steps
      if (stiff)
        ebs.set(dso::ODEStatusCode::STIFF); // Stiffness suspected
      // Copy last step
      for (int i = 0; i < neqn; i++)
        y[i] = yy(i);
      t = this->x;
      told = t;
      this->oldPermit = true;
      return ebs.set(dso::ODEStatusCode::NUMSTEPS); // Weak failure exit
    }

    // Limit step size, set weight vector and take a step
    this->h = std::copysign(
        std::min(std::abs(this->h), std::abs(tend - this->x)), this->h);
    for (int l = 0; l < neqn; l++)
      this->wt(l) = releps * std::abs(this->yy(l)) + abseps;
    int crash = this->step(this->x, this->yy(), eps);

    // Test for too small tolerances
    if (crash) {
      this->relerr = eps * releps; // Modify relative and absolute
      this->abserr = eps * abseps; // accuracy requirements
      for (int i = 0; i < neqn; i++)
        y[i] = yy(i);
      t = this->x;
      told = t;
      this->oldPermit = true;
      return ebs.set(dso::ODEStatusCode::BADACC); // Weak failure exit
    }

    ++nostep; // Count total number of steps

    // Count number of consecutive steps taken with the order of
    // the method being less or equal to four and test for stiffness
    ++kle4;
    if (this->kold > 4)
      kle4 = 0;
    if (kle4 >= 50)
      stiff = true;

  } // end step loop

  return ebs;
}
