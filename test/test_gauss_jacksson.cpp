#include "gauss_jackson.hpp"
#include <cmath>
#include <cstdio>

using dso::Vector3;

// r'' = -r / |r|^3
void diffeq([[maybe_unused]] double t, const Vector3 &r,
            [[maybe_unused]] const Vector3 &v, Vector3 &a) noexcept {
  const double norm = std::sqrt(r.x() * r.x() + r.y() * r.y());
  a.x() = -r.x() / std::pow(norm, 3);
  a.y() = -r.y() / std::pow(norm, 3);
}

int main() {
  const double e = 0.1;     /* constant value, eccentricity */
  const double tend = 20e0; /* target t, we want y(t) for t=tend */
  /* different number of iterations per example solution */
  const int Steps[] = {50, 100, 250, 500, 750, 1000, 1500, 2000};
  /* elements in Steps array */
  int numSteps = sizeof(Steps) / sizeof(Steps[0]);

  /* auxiliary */
  double t0 = 0e0, h;
  Vector3 r0{{1e0 - e, 0e0, 0e0}};
  Vector3 v0{{0e0, std::sqrt((1e0 + e) / (1e0 - e)), 0e0}};

  /* loop through iteration numbers ... */
  for (int j = 0; j < numSteps; j++) {

    /* step size h */
    h = tend / Steps[j];

    /**/
    dso::GaussJacksonIntegrator<4> GJ4(h,diffeq);
    GJ4.initialize(t0,r0,v0);

    Vector3 r,v;
    double t=t0;
    for (int i=1; i<=Steps[j]; i++)
        t = GJ4.step(t, r, v);
    
    printf("\nNums Steps: %5d State Vector: [%.5f, %.5f, %.5f, %.5f, %.5f, %.5f]\n", Steps[j], r.x(), r.y(), r.z(), v.x(), v.y(), v.z());
  }

  return 0;
}
