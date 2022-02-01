#include "earth_gravity.hpp"
#include <cstdio>

int main(int argc, char *argv[]) {
    if (argc!=2) {
        fprintf(stderr, "Usage: %s [GRAVITY MODEL FILE]\n", argv[0]);
        return 1;
    }

    EarthGravityModel gmod;
    printf("reading gravity model from file %s\n", argv[1]);

    if (gmod.init(argv[1])) {
      fprintf(stderr, "Failed parsing model from file %s \n", argv[1]);
      return 1;
    }

    #ifdef DEBUG
    gmod.print();
    #endif

    return 0;
}