#include "egravity.hpp"
#include "icgemio.hpp"
#include <cstdio>

int dso::parse_gravity_model(const char *model_fn, int degree, int order,
                             const dso::TwoPartDate &t,
                             dso::StokesCoeffs &harmonics
                             ) noexcept {
  dso::Icgem gfc(model_fn);

  /* parse the header ... */
  if (gfc.parse_header()) {
    fprintf(stderr,
            "[ERROR] Failed to parse icgem header for %s (traceback: %s)!\n",
            model_fn, __func__);
    return 1;
  }

  /* check degree and order validity */
  if (!(degree <= gfc.degree() && order <= degree)) {
    fprintf(stderr,
            "[ERROR]  invalid degree/order %d/%d for input gravity model %s "
            "(traceback: %s)\n",
            degree, order, model_fn, __func__);
    return 1;
  }

  /* resize HarmonicCoeffs to fit input */
  harmonics.resize(degree, order);

  /* parse data; store coefficients to harmonics */
  if (gfc.parse_data(degree, order, t, &harmonics)) {
    fprintf(stderr,
            "[ERROR] Failed to parse harmonic coefficients from file %s "
            "(traceback: %s)\n",
            model_fn, __func__);
    return 1;
  }

  /* all done */
  return 0;
}
