#include "cweb.hpp"
#include <cstdio>

int main() {
  printf("try downloading a valid url ...\n");
  int status = dso::http_get("https://hpiers.obspm.fr/iers/bul/bulb_new/bulletinb.252", "foobar.bb");
  printf("status=%d\n", status);
  
  printf("try downloading a non-valid url ...\n");
  status = dso::http_get("https://hpiers.obspm.fr/iers/bul/bulb_new/bulletinb.abc252", "foobar.bb");
  printf("status=%d\n", status);

  return 0;
}
