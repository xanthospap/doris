#include "sp3.hpp"
#include <bits/c++config.h>
#include <cstdio>

using namespace ids;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage " << argv[0] << " SP3_FILE\n";
    return 1;
  }

  Sp3c sp3(argv[1]);
  sp3.print_members();

  SatelliteId sv("L27");
  Sp3DataBlock block;

  if (!sp3.has_sv(sv)) {
    std::cout << "Satellite \'" << sv.to_string()
              << "\' not included in sp3 file.\n";
    return 0;
  }

  // let's try reading the records; note that -1 denotes EOF
  int j;
  std::size_t rec_count = 0;
  do {
    printf("reading record #%6lu", rec_count);
    j = sp3.get_next_data_block(sv, block);
    if (j > 0) {
      printf("Something went wrong ....status = %3d\n", j);
      return 1;
    } else if (j == -1) {
      printf("EOF encountered; Sp3 file read through!\n");
    }
    bool position_ok = !block.flag.is_set(Sp3Event::bad_abscent_position);
    bool velocity_ok = !block.flag.is_set(Sp3Event::bad_abscent_velocity);
    printf(" position: %1d velocity %1d\n", (int)position_ok, (int)velocity_ok);
    ++rec_count;
  } while (!j);

  return 0;
}
