#include "nrlmsise00.hpp"

using namespace dso::nrlmsise00;

int dso::Nrlmsise00::gtd7d(const InParams *in, int mass,
                           OutParams *out) noexcept {
  gtd7(in, mass, out);

  double *__restrict__ d = out->d;

#if __cplusplus >= 202002L
  if (mass == 48) [[likely]]{
#else
  if (mass == 48) {
#endif
    d[5] = 1.66e-24 * (4e0 * d[0] + 16e0 * d[1] + 28e0 * d[2] + 32e0 * d[3] +
                       40e0 * d[4] + d[6] + 14e0 * d[7] + 16e0 * d[8]);
  }

  if (in->meters())
    d[5] /= 1e3;

  return 0;
}
