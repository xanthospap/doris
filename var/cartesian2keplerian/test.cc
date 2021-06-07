#include "orbcrd.hpp"
#include <iostream>
#include <fstream>
#include <chrono>

using namespace std::chrono;

int main(int argc, char* argv[]) {
  if (argc<2) {
    std::cerr<<"Usage: "<<argv[0]<<" SV_STATE_FILE\n";
    return 1;
  }

  double x[3], v[3];
  std::ifstream fin(argv[1]);
  std::size_t count=0;
  auto start = high_resolution_clock::now();
  while (fin && ++count) {
    fin >> x[0] >> x[1] >> x[2] >> v[0] >> v[1] >> v[2];
    // printf("%+15.6f %+15.6f %+15.6f %+15.6f %+15.6f %+15.6f\n", x[0], x[1], x[2] ,v[0], v[1], v[2]);
#ifdef SCHWARZ
    cartesian2keplerian_schwarz(x, v);
#elif BEUTLER
    cartesian2keplerian_beutler(x, v);
#elif BERNESE
    cartesian2keplerian_bernese(x ,v);
#elif TAPLEY
    cartesian2keplerian_tapley(x, v);
#endif
  }
  auto stop = high_resolution_clock::now();

  auto duration = duration_cast<microseconds>(stop - start);
#ifdef SCHWARZ
  printf("Schwarz implementation: ");
#elif BEUTLER
  printf("Beutler implementation: ");
#elif BERNESE
  printf("Bernese implementation: ");
#elif TAPLEY
  printf("Tapley implementation: ");
#endif
  std::cout << "duration: " << duration.count() << " microseconds for "<< count << " lines\n";

  return 0;
}
