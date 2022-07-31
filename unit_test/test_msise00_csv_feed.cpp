#include "atmosphere.hpp"
#include <cstdio>

// dso::nrlmsise00::InParams<dso::nrlmsise00::NONE> Results[];


int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr,
            "ERROR. Must provide CelesTrak Space-Weather CSV data file\n");
    return 1;
  }

  dso::datetime<dso::seconds> t(dso::year(2022), dso::month(7),
                                dso::day_of_month(25), dso::seconds(0));

  dso::nrlmsise00::InParams<dso::nrlmsise00::detail::FluxDataFeedType::ST_CSV_SW>
      dataHunter(argv[1], t.mjd(), 0e0);

  double sec_per_step = 60 * 60; // aka, 1 hour
  for (int step = 0; step < 5; step++) {
    printf("--------- Hours in day: %.3f ------------------\n", step * sec_per_step / (60e0*60e0));
    dataHunter.update_params(t.mjd().as_underlying_type(), step * sec_per_step);
  }

  return 0;
}
