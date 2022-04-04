#include "astrodynamics.hpp"
#include "gauss_jackson.hpp"
#include "runge_kutta_nystrom.hpp"
#include "orbit_integrators.hpp"
#include <chrono>
#include <cstdio>

using namespace std::chrono;
using dso::Vector3;

// Test numerical integrators based on the plane, 2-body problem
// r'' = - r / |r|^3
// r(0) = (1-e, 0) and r'(0) = (0, [(1+e)*(1-e)]^(1/2))

void reference_values(double t, double e, Vector3 &r, Vector3 &v) noexcept {
  const double E = dso::kepler(e, t);
  const double sE = std::sin(E);
  const double cE = std::cos(E);
  r.x() = cE - e;
  r.y() = std::sqrt(1e0 - e * e) * sE;
  r.z() = 0e0;
  v.x() = -sE / (1e0 - e * cE);
  v.y() = std::sqrt(1e0 - e * e) * cE / (1e0 - e * cE);
  v.z() = 0e0;
}

// the 2nd degree ode:
void kepler_ode([[maybe_unused]] double t, const Vector3 &r,
                [[maybe_unused]] const Vector3 &v, Vector3 &a) noexcept {
  a = (-1e0 * r) / std::pow(r.norm(), 3e0);
}

int main() {
  // geometry of orbit
  const double e = 1e-1;

  // initial conditions
  Vector3 r0{{1e0 - e, 0e0, 0e0}};
  Vector3 v0{{0e0, std::sqrt((1e0 + e) / (1e0 - e)), 0e0}};
  const double t0 = 0e0;
  const double tend = 20e0;

  // num of steps for each test
  const int Steps[] = {10, 100, 500, 1000, 10000};

  Vector3 r, v, rref, vref;
  printf("%8s %8s %8s %11s %11s %11s %11s %15s\n", "Method", "Time", "Steps",
         "DeltaRx", "DeltaRy", "DeltaVx", "DeltaVy", "MicroSec");
  printf("-------------------------------------------------------------------------------------------\n");

  for (int test = 0; test < int(sizeof(Steps) / sizeof(int)); test++) {
    double h = (tend - t0) / Steps[test];

    auto start = high_resolution_clock::now();
    dso::GaussJacksonIntegrator<4> gj4(h, kepler_ode);
    gj4.initialize(t0, r0, v0);
    double t = t0;
    for (int k = 1; k < Steps[test] + 1; k++)
      t = gj4.step(t, r, v);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    reference_values(t, e, rref, vref);
    printf("%8s %8.5f %8d %.5e %.5e %.5e %.5e %15ld\n", "GJ4", t, Steps[test],
           std::abs(r.x() - rref.x()), std::abs(r.y() - rref.y()),
           std::abs(v.x() - vref.x()), std::abs(v.y() - vref.y()),
           static_cast<long>(duration.count()));
    
    start = high_resolution_clock::now();
    dso::GaussJacksonIntegrator<8> gj8(h, kepler_ode);
    gj8.initialize(t0, r0, v0);
    t = t0;
    for (int k = 1; k < Steps[test] + 1; k++)
      t = gj8.step(t, r, v);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    reference_values(t, e, rref, vref);
    printf("%8s %8.5f %8d %.5e %.5e %.5e %.5e %15ld\n", "GJ8", t, Steps[test],
           std::abs(r.x() - rref.x()), std::abs(r.y() - rref.y()),
           std::abs(v.x() - vref.x()), std::abs(v.y() - vref.y()),
           static_cast<long>(duration.count()));
    
    start = high_resolution_clock::now();
    dso::RungeKuttaNystromIntegrator<4,3,dso::orbit_integrators::StepSizeControl::Off> rkn43o(h, kepler_ode);
    t = t0;
    for (int k = 1; k < Steps[test] + 1; k++)
      t = rkn43o.step(t, r, v);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    reference_values(t, e, rref, vref);
    printf("%8s %8.5f %8d %.5e %.5e %.5e %.5e %15ld\n", "RKN43O", t, Steps[test],
           std::abs(r.x() - rref.x()), std::abs(r.y() - rref.y()),
           std::abs(v.x() - vref.x()), std::abs(v.y() - vref.y()),
           static_cast<long>(duration.count()));
    
    start = high_resolution_clock::now();
    dso::RungeKuttaNystromIntegrator<6,4,dso::orbit_integrators::StepSizeControl::Off> rkn64o(h, kepler_ode);
    t = t0;
    for (int k = 1; k < Steps[test] + 1; k++)
      t = rkn64o.step(t, r, v);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    reference_values(t, e, rref, vref);
    printf("%8s %8.5f %8d %.5e %.5e %.5e %.5e %15ld\n", "RKN64O", t, Steps[test],
           std::abs(r.x() - rref.x()), std::abs(r.y() - rref.y()),
           std::abs(v.x() - vref.x()), std::abs(v.y() - vref.y()),
           static_cast<long>(duration.count()));
    
    /*start = high_resolution_clock::now();
    dso::RungeKuttaNystromIntegrator<4,3,dso::orbit_integrators::StepSizeControl::Auto> rkn43a(h, kepler_ode, 1e-6);
    t = t0;
    // for (int k = 1; k < Steps[test] + 1; k++)
    while (t<tend)
      t = rkn43a.step(t, r, v);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    reference_values(t, e, rref, vref);
    printf("%8s %8.5f %8d %.5e %.5e %.5e %.5e %15ld\n", "RKN43A", t, Steps[test],
           std::abs(r.x() - rref.x()), std::abs(r.y() - rref.y()),
           std::abs(v.x() - vref.x()), std::abs(v.y() - vref.y()),
           static_cast<long>(duration.count()));*/
  }

  return 0;
}