#include "sp3.hpp"
#include <bits/c++config.h>
#include <cstdio>

using namespace ids;

int main(int argc, char* argv[]) {
  if (argc!=2) {
    std::cerr << "Usage "<<argv[0]<<" SP3_FILE\n";
    return 1;
  }

  Sp3c sp3(argv[1]);
  sp3.print_members();

  // let's try reading the records; note that -1 denotes EOF
  int j;
  std::size_t rec_count = 0;
  do {
    printf("reading record #%6lu\n", rec_count);
    j = sp3.get_next_data_block();
    if (j>0) {
      printf("Something wen wrong ....status = %3d\n", j);
      return 1;
    } else if (j==-1) {
      printf("EOF encountered; Sp3 file read through!\n");
    }
    ++rec_count;
  } while (!j);

  return 0;
}
