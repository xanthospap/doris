#include "sgode.hpp"
#include <limits>
#ifdef DEBUG
#include <cstdio>
#endif

/*
 * Suppose we have a system of first order DE:
 *             y'_1(t) = f_1(t, y_1, y_2, ... ,y_n)
 *             y'_2(t) = f_2(t, y_1, y_2, ... ,y_n)
 *             ........
 *             y'_n(t) = f_n(t, y_1, y_2, ... ,y_n)
 * with initial conditions           
 *             y_1(a), y_2(a), ..., y_n(a)
 * which can be (and hereafter will be) represented in vector form:
 *             y' = f(t,y(t))                                              (1)
 *               y(a) = y_0                                                (2)
 * we want a solution of the system at t=b, that is y(b).
 *
 * To solve the problem we need:
 *  1. define the equations (usually a function that evaluates them)
 *  2. initial conditions
 *  3. definition of the time interval [a,b]
 *  4. declaring what accuracy is expected and how the error is to be measured
 *
 *  If the integration/program fails, it should report the solutions at the 
 *  place that it failed and why it failed. DE is intended to be used for 
 *  problems in which only the solutions at the endpoint b are of interest or, 
 *  more commonly, when a table of the solutions at a sequence of output 
 *  points is desired.
 *
 *  All essential information is passed to and from DE with only eight 
 *  parameters:
 *           SUBROUTINE DE(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG)
 *  where:
 *    1. NEQN represents the number of equations n in (1); in our case
 *       this is the instance's neqn variable
 *    2. F is the name of a subroutine defining the differential equations; in
 *       our case, this correspond to the instance's f variable (of type 
 *       ODEfun). Note that in our case the instance params variable (a pointer 
 *       of type dso::IntegrationParameters is also used to assist evaluation 
 *       of f). A call to f, should evaluate the right-hand side of (1)
 *    3. T and TOUT (or in out case t and tout) represent the time interval 
 *       [a,b] (b<a can hold with no problem).
 *  The program can change the values of the calling instance's variables:
 *    * told,
 *    * relerr,abserr and
 *    * iflag
 *
 * The code attempts at each internal step to control each component of the 
 * local error vector so that
 *                 | LocalError_k | <= relerr * |y(k)| + abserr
 * This is a mixed relative-absolute error criterion that includes, as special 
 * cases, pure absolute error (RELERR = 0) and pure relative error 
 * (ABSERR = 0).
 * We suggest the following as rules of thumb for choosing a criterion. If 
 * the solution changes a great deal in magnitude during the integration, and 
 * you wish to see this change, use relative error. If the solution does not 
 * vary much, or if you are not interested in it when it is small, use 
 * absolute error. A mixed criterion is probably the best and safest choice. 
 * For solutions large in magnitude it is essentially relative error and for 
 * solutions small in magnitude it is essentially absolute error. Thus it is 
 * always a reasonable choice and avoids the troubles of the pure criteria.
 *
 * A full explanation of the IFLAG (status of the function call), is given in 
 * the respective header file.
 *
 * There is a potential source of difficulty which is dealt with by a system 
 * of negative values of IFLAG. To use Adams methods efficiently, we must let 
 * the code step along using the largest steps that will yield the requested 
 * accuracy everywhere. Between mesh points the answers are determined by 
 * interpolation very cheaply. If the step sizes chosen by the code are 
 * smaller than the interval of integration, there is no difficulty. But if 
 * the user wants a lot of output, the interval of integration may be 
 * substantially smaller than the step size the code could use. To handle this 
 * possibility efficiently the code must be allowed to integrate past TOUT 
 * internally and then compute the answers at TOUT by interpolation. The 
 * farther it is permitted to go, the more robust the code is with respect to 
 * a lot of output. For nearly all problems it is permissible to integrate 
 * past TOUT, though obviously some limit must be imposed. Our experience has 
 * been that users rarely ask for output at a spacing less than one tenth the 
 * natural step size, though it is quite common for them to want it more often 
 * than at mesh points alone. To deal with this DE permits the integration to 
 * go past TOUT internally, though never past T + 10*(TOUT — T). For nearly 
 * all problems this design means that the cost is insensitive to the number 
 * and placement of the TOUT and the user need give the matter no thought at 
 * all. For some problems it is not permissible for DE to go past TOUT 
 * internally, e.g., the solution or its derivative is not defined beyond 
 * TOUT, or there is a discontinuity there. To warn the code of this, and so 
 * to have DE stop internally at TOUT, set IFLAG negative. Thus to initialize 
 * the code and have it stop at TOUT internally, use IFLAG = — 1. When 
 * continuing an integration which is to be stopped internally at TOUT, use 
 * IFLAG = — 2. Since the code is designed to return to the user in the event 
 * of trouble with a diagnostic in IFLAG, and be set to continue with a simple 
 * call again, there are negative values of IFLAG corresponding to the returns 
 * of IFLAG = 3, 4, 5. Thus if DE does not reach TOUT and returns with, for 
 * example, IFLAG = — 3, then on calling DE again it will try once more to 
 * reach TOUT and stop internally at that point.
 *
 * The only legitimate situation in which one would want to stop at TOUT and 
 * then to continue the integration is in dealing with a discontinuity, or the 
 * like. Such a situation ought to be followed by a restart.
 */

/// According to Shampine & Gordon:
/// "The only machine dependent constants are TWOU and FOURU, which
/// represent, respectively, two times and four times the machine's unit 
/// roundoff errror U. This latter quantity is defined as being the smallest 
/// positive number for which 1+U > 1."
/// So, i guess in our case we define:
/// U = umach = std::numeric_limits<double>::epsilon();
constexpr const double umach = std::numeric_limits<double>::epsilon();
constexpr const double fouru = 4e0 * umach;

/*
dso::SGOde::IFLAG dso::SGOde::de_start() noexcept {
  if (neqn < 1) 
    return SGOde::IFLAG::INVALID_INPUT;
  
  double eps = std::max(relerr, abserr);
  if ((relerr < 0e0 || abserr < 0e0) || eps < 0e0)
    return SGOde::IFLAG::INVALID_INPUT;


}
*/

int dso::SGOde::de(double &t, double tout, const Eigen::VectorXd &y0,
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
  int crash = false;

  // test for improper parameters
  double eps = std::max(relerr, abserr);
  if (t == tout || (relerr < 0e0 || abserr < 0e0) || eps < 0e0 ||
      !iflag || (t != told && std::abs(iflag) != 1)) {
    iflag = 6;
#ifdef DEBUG
    int error_flag = -1;
    if (neqn < 1)
      error_flag = 1;
    else if (t == tout)
      error_flag = 2;
    else if (relerr < 0e0 || abserr < 0e0)
      error_flag = 3;
    else if (eps < 0e0)
      error_flag = 4;
    else if (!iflag)
      error_flag = 5;
    else if (t != told && std::abs(iflag) != 1)
      error_flag = 6;
    fprintf(stderr, "ERROR Invalid parameters to %s. error-flag: %d\n",
            __func__, error_flag);
#endif
    return 1;
  }

  const int isn = std::copysign(1, iflag);
  iflag = std::abs(iflag);

  if (iflag < 0 || iflag > 5)
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
      intrp(tout, yout);
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
      f(x, yy(), yp(), *params); // derivate at yp()
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
