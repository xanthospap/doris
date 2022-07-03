#include "sgode.hpp"
#include <limits>
#ifdef DEBUG
#include <cstdio>
#endif

constexpr const double umach = std::numeric_limits<double>::epsilon();
constexpr const double twou = 2e0 * umach;
constexpr const double fouru = 4e0 * umach;

int dso::SGOde::de(double& t, double tout, const Eigen::VectorXd &y0,
                   Eigen::VectorXd &yout) noexcept {
  // ***********************************************************************
  // *  the only machine dependent constant is based on the machine unit   *
  // *  roundoff error  u  which is the smallest positive number such that *
  // *  1.0+u .gt. 1.0 .  u  must be calculated and  fouru=4.0*u  inserted *
  // *  in the following data statement before using  de .  the routine    *
  // *  machin  calculates  u .  fouru  and  twou=2.0*u  must also be      *
  // *  inserted in subroutine  step  before calling  de .                 *
  // *    data fouru/.888d-15/                                           ***
  // ***********************************************************************
  //
  //    the constant  maxnum  is the maximum number of steps allowed in one
  //    call to  de .  the user may change this limit by altering the
  //    following statement
  constexpr const int maxnum = 500;

  int crash = false;

  // test for improper parameters
  double eps = std::max(relerr, abserr);
  if (neqn < 1 || t == tout || (relerr < 0e0 || abserr < 0e0) || eps < 0e0 ||
      !iflag || (t != told && std::abs(iflag) != 1)) {
    iflag = 6;
#ifdef DEBUG
    int error_flag = -1;
    if (neqn < 1)
      error_flag = 1;
    else if (t==tout)
      error_flag = 2;
    else if (relerr < 0e0 || abserr < 0e0)
      error_flag = 3;
    else if (eps < 0e0)
      error_flag = 4;
    else if (!iflag)
      error_flag = 5;
    else if (t!=told && std::abs(iflag) != 1)
      error_flag = 6;
    fprintf(stderr, "ERROR Invalid parameters to %s. error-flag: %d\n", __func__, error_flag);
#endif
    return 1;
  }

  const int isn = std::copysign(1, iflag);
  iflag = std::abs(iflag);

  if (iflag<0 || iflag > 5)
    return 1;

  // on each call set interval of integration and counter for number of
  // steps.  adjust input error tolerances to define weight vector for
  // subroutine  step
  // -- break point 20: --
  const double del = tout - t;
  const double absdel = std::abs(del);
  double tend = t + 10e0 * del;
  if (isn < 0)
    tend = tout;

  int nostep = 0;
  int kle4 = 0;
  int stiff = false;
  const double releps = relerr / eps;
  const double abseps = abserr / eps;

  if (delsgn * del <= 0e0 || iflag == 1) {
    // on start and restart also set work variables x and yy(*), store the
    // direction of integration and initialize the step size
    x = t;
    yy() = y0;
    // -- break point 30: --
    start = true;
    // -- break point 40: --
    delsgn = std::copysign(1e0, del);
    h = std::copysign(std::max(std::abs(tout - x), fouru * std::abs(x)),
                      tout - x);
  }

  while (true) {
    if (std::abs(x - t) >= absdel) {
      // if already past output point, interpolate and return
      // -- break point 50: --
      intrp(tout, yout, ypout());
      iflag = 2;
      t = tout;
      told = t;
      isnold = isn;
      return 0;
    }

    // if cannot go past output point and sufficiently close,
    // extrapolate and return
    if (!(isn > 0 || std::abs(tout - x) >= fouru * std::abs(x))) {
      // -- break point 60: --
      h = tout - x;
      f(x, yy(), yp(), params); // derivate at yp()
      yout = yy() + h * yp();
      iflag = 2;
      t = tout;
      told = t;
      isnold = isn;
      return 0;
    }

    // test too many steps
    if (nostep >= maxnum) {
      // -- break point 80: --
      iflag = isn * 4;
      if (stiff)
        iflag = isn * 5;
      // -- break point 90: --
      yout = yy();
      t = x;
      told = t;
      isnold = 1;
      return 1;
    }

    // limit step size, set weight vector and take a step
    // -- break point 100: --
    h = std::copysign(std::min(std::abs(h), std::abs(tend - x)), h);
    wt() = (releps * yy().cwiseAbs()).array() + abseps;
    this->step(eps, crash);

    // test for tolerances too small
    if (crash) {
      iflag = isn * 3;
      relerr = eps * releps;
      abserr = eps * abseps;
      yout = yy();
      t = x;
      told = t;
      isnold = 1;
      return 1;
    }

    // augment counter on number of steps and test for stiffness
    ++nostep;
    ++kle4;
    if (kold > 4)
      kle4 = 0;
    if (kle4 >= 50)
      stiff = true;
  }

  return 100;
}
