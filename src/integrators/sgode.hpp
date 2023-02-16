#ifndef __DSO_SGODE_ODE_HPP__
#define __DSO_SGODE_ODE_HPP__

#include "odefun.hpp"
#include "orbit_integration.hpp"

namespace dso {

class SGOde {
  /// The constant maxnum is the maximum number of steps allowed in one
  /// call to sgode_de. The user may change this limit by altering the
  /// following statement
  static constexpr const int maxnum = 500;

  
  enum class IFLAG : char {
    /// IFLAG = 2—integration successful, T is set to TOUT and Y to the
    /// solution at TOUT.
    /// All parameters in the call list are set for continuing the integration 
    /// if the user wishes to. All he has to do is define a new value TOUT and 
    /// call DE again.
    i2,
    /// IFLAG = 3—error tolerances RELERR and ABSERR are too small for the 
    /// machine being used. T is set to the point closest to TOUT reached 
    /// during the integration and Y to the solution at that point. RELERR and 
    /// ABSERR are set to larger, acceptable values. To continue with the 
    /// larger tolerances, just call DE again.
    /// The word length of the machine and the error criterion impose 
    /// limitations on the accuracy that can be obtained. Requests for too 
    /// much accuracy are detected without additional function evaluations 
    /// (calls to the subroutine F) and suitable tolerances are communicated 
    /// to the user. Because of the way control is returned, nothing is lost. 
    /// The user gets the last solution computed by the code which meets the 
    /// requested error tolerances and, if he wishes to continue with the 
    /// larger tolerances, he simply calls DE again. To integrate with the 
    /// maximum accuracy possible, simply specify tolerances that are known to 
    /// be too small and let the code increase them to an acceptable level. 
    /// The code will increase the tolerances if it is necessary, but will not 
    /// decrease them. The tolerances may be altered by the user at each call 
    /// without re-initializing.
    i3,
    /// IFLAG = 4— more than MAXNUM steps are required to reach TOUT. T is set 
    /// to the point closest to TOUT reached during the integration, and Y to 
    /// the answer at that point. To continue, just call DE again.
    i4,
    /// IFLAG = 5—more than MAXNUM steps needed to reach TOUT and the 
    /// equations appear to be stiff. T is set to the point closest to TOUT 
    /// reached during the integration, and Y to the answer at that point. A 
    /// code for stiff equations should be used but one can (usually) get 
    /// accurate results with DE if he is prepared to stand the cost. To 
    /// continue, the user has only to call DE again.
    i5,
    /// IFLAG = 6— integration is not begun because the input parameters are 
    /// invalid. The user must correct them and call DE again.
    INVALID_INPUT,
  };

public:
  SGOde(ODEfun _f, int _neqn, double rerr, double aerr,
        dso::IntegrationParameters *_params = nullptr)
      : f(_f), neqn(_neqn), iflag(1), relerr(rerr), abserr(aerr),
        params(_params) {
    Phi = Eigen::MatrixXd(neqn, 16);
    ArraysNeqn = Eigen::MatrixXd(neqn, 5); // wt, p, ypout, yp, yy
    Arrays13 = new double[13 * 7]; // psi, alpha, beta, v, w, sig, g (col-major)
  }

  ~SGOde() noexcept {
    if (Arrays13)
      delete[] Arrays13;
  }

  int flag() const noexcept { return iflag; }
  int &flag() noexcept { return iflag; }

  int de(double &t, double tout, const Eigen::VectorXd &y0,
         Eigen::VectorXd &yout) noexcept;
  int step(double &eps, int &crash) noexcept;
  // ypout is stored in the member variable ypout
  int intrp(double xout, Eigen::VectorXd &yout/*,
            Eigen::Ref<Eigen::VectorXd> ypout*/) noexcept;

  Eigen::Ref<Eigen::VectorXd> wt() noexcept { return ArraysNeqn.col(0); }
  Eigen::Ref<Eigen::VectorXd> p() noexcept { return ArraysNeqn.col(1); }
  Eigen::Ref<Eigen::VectorXd> yy() noexcept { return ArraysNeqn.col(2); }
  // !! Warning !! //
  // Column 3 of ArraysNeqn should be 3, do not change this!
  // see the bug in sgode_step.cpp ~line 330
  Eigen::Ref<Eigen::VectorXd> yp() noexcept { return ArraysNeqn.col(3); }
  Eigen::Ref<Eigen::VectorXd> ypout() noexcept { return ArraysNeqn.col(4); }
  double &wt(int i) noexcept { return ArraysNeqn(i, 0); }
  double &p(int i) noexcept { return ArraysNeqn(i, 1); }
  double &yy(int i) noexcept { return ArraysNeqn(i, 2); }
  double &yp(int i) noexcept { return ArraysNeqn(i, 3); }
  double &ypout(int i) noexcept { return ArraysNeqn(i, 4); }

  /// get psi array coefficient (size=12)
  double &psi(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 12);
#endif
    return Arrays13[0 * 13 + i];
  }

  /// get alpha array coefficient (size=12)
  double &alpha(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 12);
#endif
    return Arrays13[1 * 13 + i];
  }

  /// get beta array coefficient (size=12)
  double &beta(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 12);
#endif
    return Arrays13[2 * 13 + i];
  }

  /// pointer to the first element in the v array (size=12)
  double *v() noexcept { return Arrays13 + 3 * 13; }

  /// get v array coefficient (size=12)
  double &v(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 12);
#endif
    return Arrays13[3 * 13 + i];
  }

  /// pointer to the first element in the w array (size=12)
  double *w() noexcept { return Arrays13 + 4 * 13; }

  /// get w array coefficient (size=12)
  double &w(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 12);
#endif
    return Arrays13[4 * 13 + i];
  }

  /// get sig array coefficient (size=13)
  double &sig(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 13);
#endif
    return Arrays13[5 * 13 + i];
  }

  /// get g array coefficient (size=13)
  double &g(int i) noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 13);
#endif
    return Arrays13[6 * 13 + i];
  }

public:
  // private:
  ODEfun f;
  int neqn;
  int iflag;
  int start, phase1, nornd, isnold, kold, k, ns;
  Eigen::MatrixXd Phi;        // dimension
  Eigen::MatrixXd ArraysNeqn; // dimension
  double *Arrays13;
  double h, hold;
  double x;
  double /*t,*/ told;
  double delsgn;
  double relerr, abserr;
  bool integrate_past_tout{true};
  /// May store a pointer to some king of parameters that are passed in the
  /// ODE function
  dso::IntegrationParameters *params{nullptr};
}; // SGOde

} // namespace dso
#endif
