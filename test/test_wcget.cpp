#include "cweb.hpp"
#include "datetime/dtfund.hpp"
#include "iers_bulletin.hpp"
#include <cassert>
#include <cstdio>
#include <random>

int main() {
  std::random_device rd;  // obtain a random number from hardware
  std::mt19937 gen(rd()); // seed the generator
  std::uniform_int_distribution<int> dstr(54831, 62502);

  printf("try downloading a valid url ...\n");
  int status = dso::http_get(
      "https://hpiers.obspm.fr/iers/bul/bulb_new/bulletinb.252", "foobar.ok");
  if (status==-1) {
    printf(">> note: file already exists!\n");
  }
  assert(status<1);

  printf("try downloading a non-valid url ...\n");
  status = dso::http_get(
      "https://hpiers.obspm.fr/iers/bul/bulb_new/bulletinb.abc252",
      "foobar.er");
  assert(status);

  printf("try downloading a valid url (csv version) ...\n");
  status = dso::http_get(
      "https://datacenter.iers.org/data/csv/bulletinb-268.csv", "foobar.csv");
  // assert(!status);

  for (int i = 0; i < 5; i++) {
    printf("------------------------------------------------------------------\n");

    int mjd = dstr(gen);
    status = dso::download_iers_bulletinb_for(mjd, "data");

    // transform MJD to year/month/day
    dso::modified_julian_day Mjd(mjd);
    auto ymd = Mjd.to_ymd();

    printf("Downloaded Bulletin B file for mjd=%d (aka YYYYMMDD=%d%02d%02d) ",
           mjd, ymd.__year.as_underlying_type(),
           ymd.__month.as_underlying_type(), ymd.__dom.as_underlying_type());
    if (!status)
      printf("success!\n");
    else
      printf("failed!\n");
  }

  return 0;
}
