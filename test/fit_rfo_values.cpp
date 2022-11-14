#include "doris_rinex.hpp"
#include <cstdio>

int main(int argc, char *argv[]) {
  if (argc<2) {
    fprintf(stderr, "Usage %s [RINEX 1] <RINEX 2> .... <RINEX N>\n", argv[0]);
    return 1;
  }

  std::vector<const char *> rnxv;
  for (int i=1; i<argc; i++)
    rnxv.push_back(argv[i]);

  dso::PolynomialModel<dso::datetime<dso::nanoseconds>> fit(1);

  if (dso::fit_relative_frequency_offset(rnxv, fit, true, 4)) {
    fprintf(stderr, "ERROR! Failed to fit RFO values!\n");
    return 2;
  }

  printf("RFO model:\n");
  printf("Reference MJD = %.9f\n", fit.xref.as_mjd());
  printf("Polynomial Coeffs:\n");
  for (int i=0; i<fit.order()+1; i++) printf("\tCoef[%2d] = %.9f\n", i, fit.cf[i]);
  printf("\n");

  return 0;
}
