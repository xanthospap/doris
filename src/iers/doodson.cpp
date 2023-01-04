#include <cstring>
#include <cctype>
#include <stdexcept>
#include <cassert>
#include "tides.hpp"

dso::DoodsonNumber::DoodsonNumber(const char *const str) {
  // check that the input string has size >=7 and that all but the 3rd char
  // are actually numeric values
  for (int i=0; i<7; i++) {
    if ( !str[i] || !(std::isdigit(*(unsigned char*)(str+i)) || i==3) )
      throw std::runtime_error(
          "[ERROR] Failed to resolve Doodson number from string\n");
  }

  // first three ints, for variables τ, s, h
  iar[0] = str[0] - '0'; 
  iar[1] = str[1] - '0';
  iar[2] = str[2] - '0';
  
  // 4th character should be either a '.' or a ','
  if (str[3] != '.' && str[3] != ',')
      throw std::runtime_error(
          "[ERROR] Failed to resolve Doodson number from string\n");

  // next three characters, for variables p, N', ps
  iar[3] = str[4] - '0';
  iar[4] = str[5] - '0';
  iar[5] = str[6] - '0';
  
  // all done
  return;
}

char *dso::DoodsonNumber::str(char *buf) const noexcept {
  int written = sprintf(buf, "%d%d%d.%d%d%d", iar[0], iar[1], iar[2], iar[3],
          iar[4], iar[5]);
  assert(written == 7);
  return buf;
}

/*
int dso::doodson2intarray(const char *const str, int *arr) noexcept {
  // check that the input string has size >=7 and that all but the 3rd char
  // are actually numeric values
  // printf("\tResolving Doodson from: [%.7s]\n", str);
  for (int i=0; i<7; i++) {
    if ( !str[i] || !(std::isdigit(*(unsigned char*)(str+i)) || i==3) )
      return 1;
  }

  // first three ints, for variables τ, s, h
  arr[0] = str[0] - '0'; 
  arr[1] = str[1] - '0' - 5;
  arr[2] = str[2] - '0' - 5;
  
  // 4th character should be either a '.' or a ','
  if (str[3] != '.' && str[3] != ',') return 2;

  // next three characters, for variables p, N', ps
  arr[3] = str[4] - '0' - 5;
  arr[4] = str[5] - '0' - 5;
  arr[5] = str[6] - '0' - 5;
  //
  //printf("\t                      : [%d%d%d.%d%d%d]\n",arr[0],arr[1],arr[2],arr[3],arr[4],arr[5]);

  // all done
  return 0;
}
*/
