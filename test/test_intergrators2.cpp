#include "gauss_jackson.hpp"
#include "iers2010/matvec.hpp"
#include <cmath>
#include <cstdio>

using dso::Vector3;

void ode21(double t, const Vector3 &x, [[maybe_unused]] const Vector3 &v,
           Vector3 &acc) noexcept {
  acc.x() = -x.x() + std::cos(4e0 * t);
}

double analytical(double t, double x0, double v0) {
  return v0 * std::sin(t) + (x0 + (1e0 / 15e0)) * std::cos(t) -
         (1e0 / 15e0) * std::cos(4e0 * t);
}

void ode2to1(double t, const Vector3 &state, Vector3 &state_deriv) {
  state_deriv.x() = state.y();
  state_deriv.y() = -state.x() + std::cos(4e0 * t);
}
double rungeKutta4_step(double t, double h, const Vector3 &x0, Vector3 &xt) {
  Vector3 k1, k2, k3, k4, x;
  ode2to1(t, x0, k1);
  x = x0 + (h / 2e0) * k1;
  ode2to1(t, x, k2);
  x = x0 + (h / 2e0) * k2;
  ode2to1(t, x, k3);
  x = x0 + h * k3;
  ode2to1(t, x, k4);
  xt = x0 + (h / 6e0) * (k1 + 2e0 * k2 + 2e0 * k3 + k4);
  return t + h;
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
           steps[i], t, tend / steps[i], anls, r.x(), std::abs(anls - r.x()));

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
