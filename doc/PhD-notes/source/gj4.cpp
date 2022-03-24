#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>

/*
 * named type alias to simplify declarations
 * normally, the parameters are:
 * t aka the independent variable
 * position vector at t
 * velocity vector at t
 * acceleration vector [out]
 */
using ode = void(double, const double*, const double*, double*);

constexpr const double SCdeltaCoeffs[] = {
  1e0, -1e0, 1e0 / 12e0, 0e0, -1e0 / 240e0, -1e0 / 240e0, -221e0 / 60480e0, -19e0 / 60480e0, -9829e0 / 3628800e0
};

constexpr const double AMgammaCoeffs[] = {
  1e0, -0.5e0, -1e0 / 12e0, -1e0 / 24e0, -19e0 / 720e0, -3e0 / 160e0, -863e0 / 60480e0, -275e0 / 24192e0, -33953e0 / 3628800e0
};

/*
 * Vector/scalar multiplication; x, y are vectors of size 3 and s is scalar.
 * Store x + s * y in the output array out
 */
double *vecmul(const double* x, const double* y, double s, int sz, double* out)
{
  for (int i = 0; i < sz; i++)
    out[i] = x[i] + s * y[i];
  return out;
}

/* dead simple RK4 integrator */
void rk4_integrate(double &t, double h, const double *pos, const double *vel, double *acc, ode *deriv) {
  double k1[3], k2[3], k3[3], k4[3], rn[3], vn[3], xt[3], yt[3], r[3], v[3];
  std::memcpy(r, pos, sizeof(double)*3);
  std::memcpy(v, vel, sizeof(double)*3);
  
  /* get k1 */
  deriv(t,r,v,k1);
  /* get k2 */
  vecmul(rn, v, h / 2e0, 3, xt);
  vecmul(vn, k1, h / 2e0, 3, yt);
  deriv(t+h/2e0, xt, yt, k2);
  /* get k3 */
  vecmul(xn, k2, h / 2e0, 3, xt);
  vecmul(yn, k2, h / 2e0, 3, yt);
  deriv(t+h/2e0, xt, yt, k3);
  /* get k4 */
  vecmul(xn, k3, h / 2e0, 3, xt);
  vecmul(yn, k3, h / 2e0, 3, yt);
  deriv(t+h/2e0, xt, yt, k3);


}

/* Order must be in the rage [0,7] */
template <int Order>
class GaussJacksonIntegrator {
  public:
  GaussJacksonIntegrator(double t0, double step_size, const double* x0, const double* y0, ode* ode_function)
      : t(t0)
      , h(step_size)
      , f(ode_function)
  {
    std::memcpy(state[0], x0, sizeof(double) * 3);
    std::memcpy(state[3], y0, sizeof(double) * 3);
  }

  /*
   * Initialization of backwards differences from initial conditions. We will
   * leave this task to an RK integrator
   */
  void initialize()
  {
    /* first compute state deriv. at t0 */
    ode(t, state[0], state[3], bckdf_t[0]);
    for (int i=1;i<Order; i++) {

  }

  private:
  double t;               /* current time */
  double h;               /* step size */
  double state[6];        /* the state vector */
  ode* f;                 /* function pointer to the ODE */
  double bckdf_t[Order];  /* Backward differences of acceleration at t */
  double bckdf_ht[Order]; /* Backward differences of acceleration at t+h */
};

int main()
{
  const double e = 0.1;     /* constant value, eccentricity */
  const double tend = 20e0; /* target t, we want y(t) for t=tend */
  /* different number of iterations per example solution */
  const int Steps[] = { 50, 100, 250, 500, 750, 1000, 1500, 2000 };
  /* elements in Steps array */
  int numSteps = sizeof(Steps) / sizeof(Steps[0]);

  /* auxiliary */
  double t, h;
  double y[4], y0[] = { 0e0, 0e0, 0e0, 0e0 };

  /* loop through iteration numbers ... */
  for (int j = 0; j < numSteps; j++) {

    /* step size h */
    h = tend / Steps[j];

    /* initial values (state vector at t=t0) */
    y0[0] = 1e0 - e;
    y0[1] = 0e0;
    y0[2] = 0e0;
    y0[3] = std::sqrt((1e0 + e) / (1e0 - e));
  }
