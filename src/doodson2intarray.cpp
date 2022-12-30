#include <cstring>
#include <cctype>
#include "tides.hpp"

int dso::doodson2intarray(const char *const str, int *arr) noexcept {
  // check that the input string has size >=7 and that all but the 3rd char
  // are actually numeric values
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

  // all done
  return 0;
}
