#include "planetpos.hpp"
#include <cstring>

namespace {
constexpr const int FILLEN = 128;
constexpr const int TYPLEN = 32;
constexpr const int SRCLEN = 128;
} // namespace

int dso::cspice::load_if_unloaded_lsk(const char *lsk_kernel) noexcept {
  SpiceInt count, which, handle;
  SpiceChar file[FILLEN];
  SpiceChar filtyp[TYPLEN];
  SpiceChar source[SRCLEN];
  SpiceBoolean found;

  // get number of all ascii kernels loaded ... (to find if we already have
  // the LSK kernel loaded)
  ktotal_c("ALL", &count);

  if (count) {
    for (which = 0; which < count; which++) {
      kdata_c(which, "ALL", FILLEN, TYPLEN, SRCLEN, file, filtyp, source,
              &handle, &found);
      printf(">> Loaded kernel: %s [%s]\n", source, __func__);
      if (!std::strncmp(file, lsk_kernel, std::strlen(lsk_kernel)))
        return 0;
    }
  }

  printf(">> Loading kernel %s [%s]\n", lsk_kernel, __func__);
  furnsh_c(lsk_kernel);
  return 0;
}

int dso::cspice::load_if_unloaded_spk(const char *spk_kernel) noexcept {
  SpiceInt count, which, handle;
  SpiceChar file[FILLEN];
  SpiceChar filtyp[TYPLEN];
  SpiceChar source[SRCLEN];
  SpiceBoolean found;

  // get number of spk kernels loaded ...
  ktotal_c("ALL", &count);

  if (count) {
    for (which = 0; which < count; which++) {
      kdata_c(which, "ALL", FILLEN, TYPLEN, SRCLEN, file, filtyp, source,
              &handle, &found);
      printf(">> Loaded kernel: %s [%s]\n", source, __func__);
      if (!std::strncmp(file, spk_kernel, std::strlen(spk_kernel)))
        return 0;
    }
  }

  printf(">> Loading kernel %s [%s]\n", spk_kernel, __func__);
  furnsh_c(spk_kernel);
  return 0;
}

int dso::cspice::load_if_unloaded_pck(const char *pck_kernel) noexcept {
  SpiceInt count, which, handle;
  SpiceChar file[FILLEN];
  SpiceChar filtyp[TYPLEN];
  SpiceChar source[SRCLEN];
  SpiceBoolean found;

  // get number of spk kernels loaded ...
  ktotal_c("ALL", &count);

  if (count) {
    for (which = 0; which < count; which++) {
      kdata_c(which, "ALL", FILLEN, TYPLEN, SRCLEN, file, filtyp, source,
              &handle, &found);
      printf(">> Loaded kernel: %s [%s]\n", source, __func__);
      if (!std::strncmp(file, pck_kernel, std::strlen(pck_kernel)))
        return 0;
    }
  }

  printf(">> Loading kernel %s [%s]\n", pck_kernel, __func__);
  furnsh_c(pck_kernel);
  return 0;
}

int dso::get_sun_moon_GM(const char *pck_kernel, double &GMSun,
                         double &GMMoon) noexcept {
  if (pck_kernel)
    dso::cspice::load_if_unloaded_pck(pck_kernel);

  int n;
  bodvrd_c("SUN", "GM", 1, &n, &GMSun);   // [km^3/ sec^2]
  bodvrd_c("MOON", "GM", 1, &n, &GMMoon); // [km^3/ sec^2]

  return 0;
}
