#include "sp3.hpp"

using namespace ids;

int main(int argc, char* argv[]) {
  if (argc!=2) {
    std::cerr << "Usage "<<argv[0]<<" SP3_FILE\n";
    return 1;
  }

  Sp3c sp3(argv[1]);
  sp3.print_members();

  return 0;
}
