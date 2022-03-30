#include "gauss_jackson.hpp"
#include "iers2010/matvec.hpp"
#include <cmath>
#include <cstdio>

using dso::Vector3;

void ode21(double t, const Vector3 &x, const Vector3 &v,
           Vector3 &acc) noexcept {
  acc.x() = v.x() / 4e0 - x.x() + std::cos(4e0*t);
}

void ode2to1(double t, const Vector3 &state0, Vector3 &state_deriv) {
  state_deriv.x() = state0.y();
  state_deriv.y() = state0.y() / 4e0 - state0.x() + std::cos(4e0*t);
}
double rungeKutta4_step(double t, double h, const Vector3 &state0,
                        Vector3 &statet) {
  Vector3 k1, k2, k3, k4, x;
  ode2to1(t, state0, k1);
  x = state0 + (h / 2e0) * k1;
  ode2to1(t, x, k2);
  x = state0 + (h / 2e0) * k2;
  ode2to1(t, x, k3);
  x = state0 + h * k3;
  ode2to1(t, x, k4);
  statet = state0 + (h / 6e0) * (k1 + 2e0 * k2 + 2e0 * k3 + k4);
  return t + h;
}

double analytical(double t, double x0, double v0) {
  const double c378 = std::sin(3e0 * std::sqrt(7e0) * t / 8e0);
  const double c15226 = 15e0 / 226e0;
  const double c837 = 8e0 / (3e0 * std::sqrt(7e0));
  double xx =
      (c837 * (v0 + (2e0 / 113e0) - (1e0 / 8e0) * (x0 + c15226)) * c378 +
       (x0 + c15226) * c378) *
          std::exp(t / 8e0) -
      std::sin(4e0 * t) / 226e0 - std::cos(4e0 * t) * c15226;
  return xx;
}

int main() {
  const double t0 = 0;
  const double tend = 25e0;
  const Vector3 r0{{0e0, 0e0, 0e0}};
  const Vector3 v0{{0e0, 0e0, 0e0}};
  const int steps[] = {2,   5,   10,   20,   50,   100,
                       200, 500, 1000, 2000, 5000, 10000};
  const int num_steps = sizeof(steps) / sizeof(int);

  for (int i = 0; i < num_steps; i++) {
    Vector3 r, v;
    
    dso::GaussJacksonIntegrator<4> gj4(tend / steps[i], ode21);
    gj4.initialize(t0, r, v);
    double t = t0;
    for (int k = 1; k < steps[i] + 1; k++)
      t = gj4.step(t, r, v);

    double anls = analytical(t, r0.x(), v0.x());
    printf("Steps: %10d, time: %6.3f step: %6.3f Analytical: %10.5f GJ4 %10.5f "
           "Diff: %.5e\n",
           steps[i], t, tend / steps[i], anls, r.x(),
           std::abs(anls - r.x()));

    t = t0;
    Vector3 state0 = r0, state;
    for (int k = 1; k < steps[i] + 1; k++) {
      t = rungeKutta4_step(t, tend / steps[i], state0, state);
      state0 = state;
    }
    anls = analytical(t, r0.x(), v0.x());
    printf("Steps: %10d, time: %6.3f step: %6.3f Analytical: %10.5f RK4 %10.5f "
           "Diff: %.5e\n",
           steps[i], t, tend / steps[i], anls, state.x(),
           std::abs(anls - state.x()));
  }

  return 0;
}
