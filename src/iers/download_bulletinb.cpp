#include "cweb.hpp"
#include "iers_bulletin.hpp"
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <system_error>

namespace fs = std::filesystem;

int bulletinb_download_and_check(long mjd, const char *remote,
                                 /*const char *local*/ const fs::path &local,
                                 bool &deletefn) noexcept {
  int status;
  if ((status = dso::http_get(remote, local.c_str())) > 0) {
    fprintf(
        stderr,
        "ERROR. Failed to download IERS Bulletin B file %s (traceback: %s)\n",
        remote, __func__);
    return 2;
  }

  deletefn = !status;

  // let's see if we guessed right and the date is within the file downloaded
  // initialize a Bulletin B instance
  dso::IersBulletinB bb(local.c_str());

  // section 1 block, to hold file records
  dso::IersBulletinB_Section1Block block;

  // consider the 'final' values only, not preliminery
  return bb.get_section1_at(mjd, block, false);
}

int dso::download_iers_bulletinb_for(long mjd, char *localfn,
                                     const char *dir) noexcept {
  // construct the directory to download files to
  fs::path local = (dir) ? dir : fs::current_path();
  if (!fs::is_directory(local)) {
    fprintf(stderr, "ERROR. %s is not a valid directory (traceback: %s)\n",
            local.c_str(), __func__);
    return -1;
  }

  char remote[64], fnbuf[64];
  int bnumber = 253;
  long bmjd = 54831;
  if (mjd < bmjd) {
    fprintf(
        stderr,
        "ERROR. Failed to download IERS Bulletin B file for MJD %ld; given "
        "date is prior to first Bulletin file for MJD %ld (traceback: %s)\n",
        mjd, bmjd, __func__);
    return 1;
  }

  // let's make a wild guess ... each bulletin B has values for 31 days;
  int guess = (mjd - bmjd) / 31 + bnumber;
  int initial_guess = guess;

  // auxiliary fs
  bool deletetmp = false;
  std::error_code ec;

  // make the filename download and check the file
  std::sprintf(remote, "https://hpiers.obspm.fr/iers/bul/bulb_new/bulletinb.%d",
               guess);
  std::sprintf(fnbuf, "bulletinb.%d", guess);
  local /= fs::path(fnbuf);
  int status = bulletinb_download_and_check(mjd, remote, local, deletetmp);

  // error parsing file
  if (status > 0) {
    fprintf(stderr,
            "ERROR. Failed to parse IERS Bulletin B file %s (traceback: %s)\n",
            local.c_str(), __func__);
    return status;
  }

  // do we have the correct file ?
  if (!status) {
    // if needed set the filename of the downloaded file
    if (localfn)
      std::strcpy(localfn, local.string().c_str());
    return 0;
  }

  int max_tries = 5;

  // if the file is ahead of the given date ...
  if (status == dso::bulletin_details::FILE_IS_AHEAD) {
    while ((status == dso::bulletin_details::FILE_IS_AHEAD) &&
           (initial_guess - guess < max_tries)) {
      if (deletetmp)
        fs::remove(local, ec);
      --guess;
      std::sprintf(remote,
                   "https://hpiers.obspm.fr/iers/bul/bulb_new/bulletinb.%d",
                   guess);
      std::sprintf(fnbuf, "bulletinb.%d", guess);
      local.replace_filename(fs::path(fnbuf));
      status = bulletinb_download_and_check(mjd, remote, local, deletetmp);
    }
    // if it is prior ...
  } else if (status == dso::bulletin_details::FILE_IS_PRIOR) {
    while ((status == dso::bulletin_details::FILE_IS_PRIOR) &&
           (guess - initial_guess < max_tries)) {
      if (deletetmp)
        fs::remove(local, ec);
      ++guess;
      std::sprintf(remote,
                   "https://hpiers.obspm.fr/iers/bul/bulb_new/bulletinb.%d",
                   guess);
      std::sprintf(fnbuf, "bulletinb.%d", guess);
      local.replace_filename(fs::path(fnbuf));
      status = bulletinb_download_and_check(mjd, remote, local, deletetmp);
    }
  }

#ifdef DEBUG
  printf("%s Number of guesses for the Bulletin B file %d (started with number "
         "%d, ended at %d)\n",
         __func__, std::abs(guess - initial_guess), initial_guess, guess);
#endif

  // ... so, de we have the right file?
  if (!status) {
    // if neede set the filename of the downloaded file
    if (localfn)
      std::strcpy(localfn, local.string().c_str());
    return status;
  }

  if (deletetmp)
    fs::remove(local, ec);
  fprintf(stderr,
          "ERROR. Failed to download any matching IERS Bulletin B file for MJD "
          "%ld (traceback: %s)\n",
          mjd, __func__);
  return status;
}
