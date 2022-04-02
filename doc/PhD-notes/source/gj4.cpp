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
using ode = void(double, const double *, const double *, double *);

void diffeq(double t, const double *r, const double *v, double *a) {
  const double norm = std::sqrt(r[0] * r[0] + r[1] * r[1]);
  a[0] = -r[0] / std::pow(norm, 3);
  a[1] = -r[1] / std::pow(norm, 3);
}

constexpr const double SCdeltaCoeffs[] = {1e0,
                                          -1e0,
                                          1e0 / 12e0,
                                          0e0,
                                          -1e0 / 240e0,
                                          -1e0 / 240e0,
                                          -221e0 / 60480e0,
                                          -19e0 / 60480e0,
                                          -9829e0 / 3628800e0};

constexpr const double AMgammaCoeffs[] = {1e0,
                                          -0.5e0,
                                          -1e0 / 12e0,
                                          -1e0 / 24e0,
                                          -19e0 / 720e0,
                                          -3e0 / 160e0,
                                          -863e0 / 60480e0,
                                          -275e0 / 24192e0,
                                          -33953e0 / 3628800e0};

/*
 * Vector/scalar multiplication; x, y are vectors of size 3 and s is scalar.
 * Store x + s * y in the output array out
 */
double *vecmul(const double *x, const double *y, double s, int sz,
               double *out) {
  for (int i = 0; i < sz; i++)
    out[i] = x[i] + s * y[i];
  return out;
}

/*
 * dead simple RK4 integrator. sol (at output) will contain the state vector
 * values at t+h
 */
double rk4_integrate(double t, double h, const double *pos, const double *vel,
                     double *sol, ode *deriv) {
  double p1[3], p2[3], p3[3], p4[3];
  double v1[3], v2[3], v3[3], v4[3];
  double a1[3], a2[3], a3[3], a4[3];

  p1[0] = pos[0];
  p1[1] = pos[1];
  v1[0] = vel[0];
  v1[1] = vel[1];
  printf("\tRK4::integrate initial conditions: %.5f, %.5f and %.5f, %.5f\n", p1[0], p1[1], v1[0], v1[1]);
  diffeq(t, p1, v1, a1);
  printf("\t\tk1=%.5f, %.5f\n", a1[0], a1[1]);

  vecmul(vel, a1, h / 2, 2, v2);
  vecmul(pos, v1, h / 2, 2, p2);
  printf("\t\tcalling k2 with pos=%.5f, %.5f and vel=%.5f, %.5f (h=%.5f)\n", p2[0], p2[1], v2[0], v2[1], h);
  diffeq(t + h / 2, p2, v2, a2);
  printf("\t\tk2=%.5f, %.5f\n", a2[0], a2[1]);

  vecmul(vel, a2, h / 2, 2, v3);
  vecmul(pos, v2, h / 2, 2, p3);
  diffeq(t + h / 2, p3, v3, a3);
  printf("\t\tk3=%.5f, %.5f\n", a3[0], a3[1]);

  vecmul(vel, a3, h, 2, v4);
  vecmul(pos, v3, h, 2, p4);
  diffeq(t + h, p4, v4, a4);
  printf("\t\tk4=%.5f, %.5f\n", a4[0], a4[1]);

  sol[2] = sol[5] = 0e0;
  for (int i = 0; i < 2; i++)
    sol[i] = pos[i] + (h / 6e0) * (v1[i] + 2 * v2[i] + 2 * v3[i] + v4[i]);
  for (int i = 3; i < 5; i++) {
    sol[i] = vel[i-3] + (h / 6e0) * (a1[i-3] + 2 * a2[i-3] + 2 * a3[i-3] + a4[i-3]);
    //printf("\t\tsol[%d] = %.5f + %.2f*(%.5f + 2*%.5f + 2*%.5f + %.5f)\n", i, vel[i-3], h/6, a1[i], a2[i], a3[i], a4[i]);
  }

  return t + h;
}

/* Order must be in the rage [0,7] */
template <int Order> class GaussJacksonIntegrator {
public:
  GaussJacksonIntegrator(double t0, double step_size, const double *x0,
                         const double *y0, ode *ode_function)
      : t(t0), h(step_size), f(ode_function) {
    std::memcpy(state, x0, sizeof(double) * 3);
    std::memcpy(state+3, y0, sizeof(double) * 3);
  }

private:
  double t;        /* current time */
  double h;        /* step size */
  double state[6]; /* the state vector */
  ode *f;          /* function pointer to the ODE */
  /* Order * 3(=dimension of acceleration vector) */
  double bckdf_t[Order * 3];  /* Backward differences of acceleration at t */
  double bckdf_ht[Order * 3]; /* Backward differences of acceleration at t+h */
  double S1[3], S2[3];

public:
  /*
   * Initialization of backwards differences from initial conditions. We will
   * leave this task to an RK integrator
   */
  void initialize() {
  printf("RK Initialization:\n");
    double tn = this->t;

    // Create table of accelerations at past times t-3h, t-2h, and t-h using
    // RK4 steps
    f(t, state, state+3, bckdf_t);
    printf("\tf(t=t0) = [%.5f, %.5f, %.5f]\n", bckdf_t[0], bckdf_t[1], bckdf_t[2]);
    double rk4sol[6];
    for (int i = 1; i < Order; i++) {
      tn = rk4_integrate(tn, -h, state, state+3, rk4sol, f);
      printf("\tRK4 state=[%.5f, %.5f, %.5f, %.5f, %.5f, %.5f]\n", rk4sol[0], rk4sol[1], rk4sol[2], rk4sol[3], rk4sol[4], rk4sol[5]);
      f(tn, rk4sol, rk4sol+3, bckdf_t+(i * 3));
      std::memcpy(state, rk4sol, sizeof(double)*3);
      std::memcpy(state+3, rk4sol+3, sizeof(double)*3);
    }
    
    for (int k=0; k<4; k++)
      printf("\tbckdf_tiffs[%d] = [%.5f, %.5f, %.5f]\n", k, bckdf_t[k*3], bckdf_t[k*3]+1, bckdf_t[k*3]+2);

    // compute backward differences
    for (int i = 1; i < Order; i++) {
      for (int j = Order - 1; j >= i; j--) {
        for (int k = 0; k < 3; k++) {
          bckdf_t[j * 3 + k] = bckdf_t[(j - 1) * 3 + k] - bckdf_t[j * 3 + k];
        }
      }
    }

    // Initialize backwards sums using 4th order GJ corrector
    for (int i = 0; i < 3; i++)
      S1[i] = state[3 + i] / h;
    for (int i = 1; i <= Order; i++)
      for (int j = 0; j < 3; j++)
        S1[j] -= AMgammaCoeffs[i] * bckdf_t[(i - 1) * 3 + j];
    for (int i = 0; i < 3; i++)
      S2[i] = state[i] / (h * h) - SCdeltaCoeffs[1] * S1[i];
    for (int i = 2; i <= Order + 1; i++)
      for (int j = 0; j < 3; j++)
        S2[j] -= SCdeltaCoeffs[i] * bckdf_t[(i - 2) * 3 + j];
  
  printf("\tS1=[%.5f, %.5f, %.5f]\n", S1[0], S1[1], S1[2]);
  printf("\tS2=[%.5f, %.5f, %.5f]\n", S2[0], S2[1], S2[2]);
  }

};

int main() {
  const double e = 0.1;     /* constant value, eccentricity */
  const double tend = 20e0; /* target t, we want y(t) for t=tend */
  /* different number of iterations per example solution */
  const int Steps[] = {50, 100, 250, 500, 750, 1000, 1500, 2000};
  /* elements in Steps array */
  int numSteps = sizeof(Steps) / sizeof(Steps[0]);

  /* auxiliary */
  double t0, h;
  double state0[] = {1e0-e, 0e0, 0e0, 0e0, std::sqrt((1e0+e)/(1e0-e)), 0e0};


  /* loop through iteration numbers ... */
  for (int j = 0; j < 2/*numSteps*/; j++) {

    /* step size h */
    h = tend / Steps[j];

    /**/
    GaussJacksonIntegrator<4> GJ4(t0, h, state0, state0+3, diffeq);
    GJ4.initialize();
  }

  return 0;
}
