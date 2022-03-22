#include <cmath>
#include <cstdio>
#include <cassert>

/*
 * Vector/scalar multiplication; x, y are vectors of size sz and s is scalar.
 * Store x + s * y in the output array out
 */
void vecmul(const double *x, const double *y, double s, int sz, double *out) {
  for (int i = 0; i < sz; i++)
    out[i] = x[i] + s * y[i];
}

/* 
 * Compute state derivative at t, given the state at t.
 * Note that the parameter t in this case is actually not used.
 * The ODE reads (state=[x,y,Vx,Vy] and state_derivate=[Vx,Vx,Accx,Accy]):
 * Vx = 1 * Vx
 * Vy = 1 * Vy
 * Accx = -x / (x^2+y^2)**3
 * Accy = -y / (x^2+y^2)**3
 */
void diffeq(/*double t,*/ const double *state, double *state_deriv) {
  state_deriv[0] = state[2];
  state_deriv[1] = state[3];
  const double r = std::sqrt(state[0]*state[0] + state[1]*state[1]);
  const double r3 = r * r * r;
  state_deriv[2] = -state[0] / r3;
  state_deriv[3] = -state[1] / r3;
}

/*
 * Perform the integration of the ODE diffeq using RK4 method.
 * Starting at t=t0, perform a given number of iterations with step size h.
 * xt0 is the state vector at time t0.
 * y is the result, aka the value of the state vector after the iterations
 * at time t (to be returned).
 */
double rk4(double t0, double h, int iterations, const double *xt0, double *y) {
  /* auxiliary arrays */
  double k1[4], k2[4], k3[4], k4[4], ytmp[4], state[4];

  /* note that xt0 is the initial state */
  for (int i=0; i<4; i++) state[i] = xt0[i];
  double t = t0;
  
  for (int i = 0; i < iterations; i++) {
    diffeq(/*t0,*/ state, k1); /* state deriv. stored at k1 ... */
    vecmul(state, k1, h / 2e0, 4, ytmp); /* ytmp = state + (h/2)*k1 */
    diffeq(/*t0 + h / 2,*/ ytmp, k2);
    vecmul(state, k2, h / 2e0, 4, ytmp);
    diffeq(/*t0 + h / 2,*/ ytmp, k3);
    vecmul(state, k3, h, 4, ytmp);
    diffeq(/*t0 + h,*/ ytmp, k4);

    /* update state vector at t(n+1) */
    for (int j = 0; j < 4; j++)
      state[j] += (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
    
    /* update t */
    t += h;
  }

  /* copy the values of the state vector at t to the result array y */
  for (int i=0; i<4; i++) y[i] = state[i];
  return t;
}

/*
 * Sove the ODE: y'' = y / |y|^3 using Runge-Kutte 4 method.
 * Initial conditions: y(0) = [1-e, 0] and y'(0) = [0, {(1+e)/(1-e)}^(1/2)]
 * Use different number of steps/iterations to check the accuracy.
 */
int main() {
  const double e = 0.1; /* constant value, eccentricity */
  const double tend = 20e0; /* target t, we want y(t) for t=tend */
  /* different number of iterations per example solution */
  const int Steps[] = {50, 100, 250, 500, 750, 1000, 1500, 2000};
  /* elements in Steps array */
  int numSteps = sizeof(Steps) / sizeof(Steps[0]);

  /* auxiliary */
  double t, h;
  double y[4], y0[]={0e0,0e0,0e0,0e0};

  /* loop through iteration numbers ... */
  for (int j = 0; j < numSteps; j++) {

    /* step size h */
    h = tend / Steps[j];

    /* initial values (state vector at t=t0) */
    y0[0] = 1e0 - e;
    y0[1] = 0e0;
    y0[2] = 0e0;
    y0[3] = std::sqrt((1e0 + e) / (1e0 - e));

    /* integrate via RK-4 until t=tend; get the state vector */
    t=rk4(0, h, Steps[j], y0, y);

    /* Results */
    printf("State Vector at t = %.2f: y=[", t);
    printf("%.5f, %.5f, %.5f, %.5f]\n", y[0], y[1], y[2], y[3]);
  }

  return 0;
}
