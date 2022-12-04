#include "egravity.hpp"
#include "icgemio.hpp"
#include <cstdio>

int dso::parse_gravity_model(const char *model_fn, int degree, int order,
                             const dso::datetime<dso::nanoseconds> &t,
                             dso::HarmonicCoeffs &harmonics,
                             bool denormalize) noexcept {
  dso::Icgem gfc(model_fn);

  // parse the header ...
  if (gfc.parse_header()) {
    fprintf(stderr,
            "[ERROR] Failed to parse icgem header for %s (traceback: %s)!\n",
            model_fn, __func__);
    return 1;
  }

  if (!(degree <= gfc.degree() && order <= degree)) {
    fprintf(stderr,
            "[ERROR]  invalid degree/order %d/%d for input gravity model %s "
            "(traceback: %s)\n",
            degree, order, model_fn, __func__);
    return 1;
  }

  // resize HarmonicCoeffs to fit input
  harmonics.resize(degree, order);

  // parse data; store coefficients to harmonics
  if (gfc.parse_data(degree, order, t, &harmonics)) {
    fprintf(stderr,
            "[ERROR] Failed to parse harmonic coefficients from file %s "
            "(traceback: %s)\n",
            model_fn, __func__);
    return 1;
  }

  // for (int n=0; n<=3; n++) {
  //   printf("C(%d,%d) = %.15e ", n,0,harmonics.C(n,0));
  //   for (int m=1; m<=n; m++) {
  //     printf("C(%d,%d) = %.15e, S(%d,%d) = %.15e
  //     ",n,m,harmonics.C(n,m),n,m,harmonics.S(n,m));
  //   }
  //   printf("\n");
  // }

  // if needed denormalize coefficients
  if (denormalize)
    harmonics.denormalize();

  //for (int n = 0; n <= 3; n++) {
  //  printf("C(%d,%d) = %.15e ", n, 0, harmonics.C(n, 0));
  //  for (int m = 1; m <= n; m++) {
  //    printf("C(%d,%d) = %.15e, S(%d,%d) = %.15e ", n, m, harmonics.C(n, m), n,
  //           m, harmonics.S(n, m));
  //  }
  //  printf("\n");
  //}

  // all done
  return 0;
}
