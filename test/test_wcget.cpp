#include "cweb.hpp"
#include <cstdio>
#include <cassert>

int main() {
  printf("try downloading a valid url ...\n");
  int status = dso::http_get("https://hpiers.obspm.fr/iers/bul/bulb_new/bulletinb.252", "foobar.ok");
  assert(!status);
  
  printf("try downloading a non-valid url ...\n");
  status = dso::http_get("https://hpiers.obspm.fr/iers/bul/bulb_new/bulletinb.abc252", "foobar.er");
  assert(status);

  printf("try downloading a valid url (csv version) ...\n");
  status = dso::http_get("https://datacenter.iers.org/data/csv/bulletinb-268.csv", "foobar.csv");
  assert(!status);

  return 0;
}
