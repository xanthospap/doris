#include "astrodynamics.hpp"
#include <cmath>
#ifdef ASSERT_ERROR
#include <cassert>
#endif
#ifdef DEBUG
#include <cstdio>
#endif
/*
using dso::Mat3x3;
using dso::Vector3;

int dso::kepler2state(const double *keplerian, double *state,
                      double gm) noexcept {
  // auxiliary
  const double costh = std::cos(keplerian[5]); // cos(θ)
  const double sinth = std::sin(keplerian[5]); // sin(θ)

  // Calculate position vector r_x in perifocal coordinates
  const double rmag = ((keplerian[0] * keplerian[0]) / gm) *
                      (1e0 / (1e0 + costh * keplerian[3]));
  const Vector3 rv{{rmag * costh, rmag * sinth, 0e0}};

  // Calculate the velocity vector v_x in perifocal coordinates
  const double vmag = gm / keplerian[0];
  const Vector3 vv{{-vmag * sinth, vmag * (keplerian[3] + costh), 0e0}};

  // Calculate the matrix Q_xx' of the transformation from perifocal to
  // geocentric equatorial coordinates. Q is actually the Direction Cosine
  // Matrix.
  const double si = std::sin(keplerian[1]); // trigonometric numbers for i
  const double ci = std::cos(keplerian[1]);
  const double sO = std::sin(keplerian[2]); // trigonometric numbers for Ω
  const double cO = std::cos(keplerian[2]);
  const double so = std::sin(keplerian[4]); // trigonometric numbers for ω
  const double co = std::cos(keplerian[4]);
  const Mat3x3 Q{{-sO * ci * so + cO * co, -sO * ci * co - cO * so, sO * si,
                  cO * ci * so + sO * co, cO * ci * co - sO * so, -cO * si,
                  si * so, si * co, ci}};

  // Transform r and v into the geocentric frame via Q
  Vector3 r = Q * rv;
  std::memcpy(state, r.data, sizeof(double) * 3);
  r = Q * vv;
  std::memcpy(state + 3, r.data, sizeof(double) * 3);

  return 0;
}
*/
