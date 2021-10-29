#include "doris_rinex.hpp"
#include "doris_utils.hpp"

int main(int argc, char *argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage %s [DORS RINEX1] [... [DORIS RINEX N] ...]\n", argv[0]);
    return 1;
  }

  // aim ...
  printf("Using %d RINEX files to model receiver frequency offset\n", argc);

  // allocate and fill the list of site id's
  int num_rinex = argc-1;
  char **rinex_fns = new char *[num_rinex];
  for (int i = 0; i < num_rinex; i++) {
    rinex_fns[i] = new char[std::strlen(argv[i] + 1)];
    std::strcpy(rinex_fns[i], argv[i+1]);
  }

  // try to model the RFO parameter using the values extracted from the RINEX's
    int error = ids::fit_relative_frequency_offset(rinex_fns, num_rinex, 5.0e0, 5.0e0, 10.0e0);
    printf("All done! fit function returned %d\n", error);

  // deallocate memmory
  for (int i = 0; i < num_rinex; i++)
    delete[] rinex_fns[i];
  delete[] rinex_fns;

  return error;
}