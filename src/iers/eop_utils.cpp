#include "iers_bulletin.hpp"

int dso::fill_eop_lookup_tables(const dso::datetime<dso::nanoseconds> &t, 
  char *bulletinb_fn, dso::EopInfo &eopLUT, int size, int download = 1, bool use_prelimenery_values=false) noexcept 
{
    if (download) {
    // get the corresponding bulletin b file from IERS (download)
    int error = dso::download_iers_bulletinb_for(t.mjd().as_underlying_type(),
                                                 bulletinb_fn);
    if (error) {
      fprintf(stderr, "[ERROR] Failed to download Bulletin B file for date: %ld (traceback: %s)\n",
              t.mjd().as_underlying_type(), __func__);
      return 1;
    }
  }

  // parse the Bulletin B file
  dso::IersBulletinB_Section1Block *bbblocks =
      new dso::IersBulletinB_Section1Block[size];

  dso::IersBulletinB bulb(bulletinb_fn);
  int bbblocks_size = bulb.parse_section1(bbblocks, false);
  if (bbblocks_size <= 0) {
    fprintf(stderr, "[ERROR] Failed parsing Bulletin B file %s (traceback: %s)\n", bulletinb_fn, __func__);
    delete[] bbblocks;
    return 2;
  }

  
}
