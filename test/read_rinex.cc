#include "doris_rinex.hpp"
#include <iostream>

int main(int argc, char* argv[]) {
  if (argc!=2) {
    std::cerr<<"\nUsage: "<<argv[0]<<" [RINEX_FILE]\n";
    return 1;
  }

  ids::DorisObsRinex rnx(argv[1]);
  rnx.read();

  std::cout<<"\n";
  return 0;
}
