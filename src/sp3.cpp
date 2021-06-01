#include "sp3.hpp"
#include <stdexcept>
#ifdef DEBUG
#include <iostream>
#include "ggdatetime/datetime_write.hpp"
#endif

/// Max record characters (for a navigation data block)
constexpr int MAX_RECORD_CHARS{128};

/// @details Sp3c constructor, using a filename. The constructor will
///          initialize (set) the _filename attribute and also (try to)
///          open the input stream (i.e. _istream).
///          If the file is successefuly opened, the constructor will read
///          the header and assign info.
/// @param[in] filename  The filename of the Sp3 file
ids::Sp3c::Sp3c(const char *filename)
    : __filename(filename), __istream(filename, std::ios_base::in),
      /*__satsys(SATELLITE_SYSTEM::mixed),*/ __end_of_head(0) {
  int j;
  if ((j = read_header())) {
    if (__istream.is_open())
      __istream.close();
    throw std::runtime_error("[ERROR] Failed to read Sp3 header; Error Code: " +
                             std::to_string(j));
  }
}

#ifdef DEBUG
void ids::Sp3c::print_members() const noexcept {
  std::cout << "\nfilename     :" << __filename
            << "\nVersion      :" << version__
            << "\nStart Epoch  :" << ngpt::strftime_ymd_hms(start_epoch__)
            << "\n# Epochs     :" << num_epochs__
            << "\nCoordinate S :" << crd_sys__
            << "\nOrbit Type   :" << orb_type__
            << "\nAgency       :" << agency__
            << "\nTime System  :" << time_sys__
            << "\nInterval(sec):" << interval__.to_fractional_seconds();
}
#endif
