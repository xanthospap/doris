#include "planetpos.hpp"
#include <cstring>

#define  FILLEN   128                                                           
#define  TYPLEN   32                                                            
#define  SRCLEN   128

int dso::cspice::load_if_unloaded_lsk(const char *lsk_kernel) noexcept {
  SpiceInt count, which, handle;
  SpiceChar file[FILLEN];
  SpiceChar filtyp[TYPLEN];
  SpiceChar source[SRCLEN];
  SpiceBoolean found;

  // get number of all ascii kernels loaded ... (to find if we already have
  // the LSK kernel loaded)
  ktotal_c("text", &count);

  if (count) {
    for (which = 0; which < count; which++) {
      kdata_c(which, "text", FILLEN, TYPLEN, SRCLEN, file, filtyp, source,
              &handle, &found);
      if (!std::strncmp(file, lsk_kernel, std::strlen(lsk_kernel))) return 0;
    }
  }

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
  ktotal_c("spk", &count);

  if (count) {
    for (which = 0; which < count; which++) {
      kdata_c(which, "spk", FILLEN, TYPLEN, SRCLEN, file, filtyp, source,
              &handle, &found);
      if (!std::strncmp(file, spk_kernel, std::strlen(spk_kernel))) return 0;
    }
  }

  furnsh_c(spk_kernel);
  return 0;
}
