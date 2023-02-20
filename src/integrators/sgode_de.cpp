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

dso::SGOde::IFLAG dso::SGOde::de(double &t, double tout, const Eigen::VectorXd &y0,
                   Eigen::VectorXd &yout) noexcept {
  // test for improper parameters
  double eps = std::max(relerr, abserr);
  if ((t == tout)
      || (relerr < 0e0 || abserr < 0e0) 
      || (eps < 0e0)
      || (iflag != IFLAG::RESTART && told!=t)) {
    return (iflag = IFLAG::INVALID_INPUT);
  }

  // on each call set interval of integration and counter for number of
  // steps.  adjust input error tolerances to define weight vector for
  const double del = tout - t;
  const double absdel = std::abs(del);
  const double tend =
      integrate_past_tout * (t + 10e0 * del) + (integrate_past_tout)*tout;
  // reset number of steps
  int nostep = 0;
  int kle4 = 0;
  // reset stiff signal
  int stiff = false;
  // reset error tolerances
  const double releps = relerr / eps;
  const double abseps = abserr / eps;

  // (re-)start (could be the we are reversing the integration in time)
  if ((delsgn * del <= 0e0) || (!integrate_past_tout) || (iflag == IFLAG::RESTART)) {
    // on start and restart also set work variables x and yy
    tc = t;
    yy() = y0;
    // store direction of integration
    delsgn = (-1) * (del < 0) + (1) * (del >= 0);
    // initialize step size
    h = std::copysign(std::max(std::abs(del), fouru * std::abs(t)), del);
#ifdef INTEGRATOR_CHECK
    function_calls = 0;
    printf("[INTGR] Start of integration, calls %u start %.15e to %.15e\n", function_calls, t, tout);
#endif
  }

  while (true) {
    // if already past output point, interpolate and return
    if (std::abs(tc - t) >= absdel) {
      intrp(tout, yout);
      told = t = tout;
#ifdef INTEGRATOR_CHECK
      printf("[INTGR] End of integration, calls %u\n", function_calls);
#endif
      return (iflag = IFLAG::SUCCESS);
    }

    // if cannot go past output point and sufficiently close, extrapolate and 
    // return
    if ((!integrate_past_tout) && (std::abs(tout-tc)<fouru*std::abs(tc))) {
      h = tout - tc;
      f(tc, yy(), yp(), *params);
#ifdef INTEGRATOR_CHECK
      ++function_calls;
#endif
      yout = yy() + h * yp();
      told = t = tout;
#ifdef INTEGRATOR_CHECK
      printf("[INTGR] End of integration, calls %u\n", function_calls);
#endif
      return (iflag = IFLAG::SUCCESS);
    }

    // test too many steps
    if (nostep >= maxnum) {
      yout = yy();
      told = t = tc;
      fprintf(stderr,
              "[WRNNG] Integrator says: reached max number of steps!\n");
      if (stiff)
        fprintf(stderr, "[ERROR] Probably a stiff ODE; cannot integrate!\n");
      else
        fprintf(stderr, "[WRNNG] Allowing integration past target time, and "
                        "setting stiff flag!\n");
      integrate_past_tout = true;
      return (iflag = stiff ? IFLAG::STIFF : IFLAG::MAXSTEPS_REACHED);
    }

    // limit step size, set weight vector and take a step
    h = std::copysign(std::min(std::abs(h), std::abs(tend - tc)), h);
    // WT(l) = |y(l)|*RE/EPS + AE/EPS (shampine & Gordon, p. 176)
    wt() = (releps * yy().cwiseAbs()).array() + abseps;
    int crash = this->step(eps);

    // test for tolerances too small
    if (crash) {
      relerr *= eps;
      abserr *= eps;
      yout = yy();
      told = t = tc;
      fprintf(stderr, "[WRNNG] Integrator says: tolerances too small! Allowing "
                      "integration past target time\n");
      integrate_past_tout = true;
      return (iflag = IFLAG::TOL_SMALL);
    }

    // augment counter on number of steps and test for stiffness
    ++nostep;
    kle4 = (kold > 4) ? 0 : kle4+1;
    stiff = (kle4 >= 50);
  }

  return (iflag=IFLAG::UNDEFINED);
}
