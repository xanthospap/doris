#include "var_utils.hpp"
#include <cstdio>

int main(int argc, char *argv[]) {
  if (argc!=2) {
    fprintf(stderr, "Usage %s [CSV FLUX FILE]\n", argv[0]);
    return 1;
  }

  dso::datetime<dso::seconds> tt(dso::year(2022), dso::month(1),
                                 dso::day_of_month(1), dso::seconds());

  dso::utils::celestrack::details::CelestTrackSWFlux flux_data;
  if (int error; (error=dso::get_CelesTrack_flux_data(tt, argv[1], flux_data))) {
    fprintf(stderr, "ERROR. Failed with error code: %d\n", error);
    return error;
  }

  printf("Flux data for requested date:\n");
  printf("F10.7_OBS         = %+.2f\n", flux_data.f107Obs);
  printf("F10.7_ADJ         = %+.2f\n", flux_data.f107Adj);
  printf("F10.7_OBS_CENTER81= %+.2f\n", flux_data.f107ObsC81);
  printf("F10.7_OBS_LAST81  = %+.2f\n", flux_data.f107ObsL81);
  printf("F10.7_ADJ_CENTER81= %+.2f\n", flux_data.f107AdjC81);

  return 0;
}
