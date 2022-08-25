#ifndef __DSO__DORIS_BEACON_TABLE_IN_HPP__
#define __DSO__DORIS_BEACON_TABLE_IN_HPP__

#include "datetime/datetime_read.hpp"
#include <cstring>
#include <fstream>
#include <vector>

namespace dso {

struct BeaconInformationTableEntry {
  char _4charid[4];
  double _height;
  dso::datetime<dso::nanoseconds> from{dso::datetime<dso::nanoseconds>::max()};
  dso::datetime<dso::nanoseconds> to{dso::datetime<dso::nanoseconds>::max()};
  BeaconInformationTableEntry(const char *site, double hgt) noexcept
      : _height(hgt) {std::memcpy(_4charid,site,4*sizeof(char));
  }
  BeaconInformationTableEntry(
      const char *site, double hgt,
      const dso::datetime<dso::nanoseconds> &tfrom,
      const dso::datetime<dso::nanoseconds> &tto) noexcept
      : _height(hgt), from(tfrom), to(tto) {
    std::memcpy(_4charid, site, 4 * sizeof(char));
  }
};

class BeaconInformationTable {
  char fn_[256];

public:
  BeaconInformationTable(const char *fn) noexcept { std::strcpy(fn_, fn); }
  BeaconInformationTable(const BeaconInformationTable &) = delete;
  BeaconInformationTable(BeaconInformationTable &&) = delete;
  BeaconInformationTable &operator=(const BeaconInformationTable &) = delete;
  BeaconInformationTable &operator=(BeaconInformationTable &&) = delete;
  ~BeaconInformationTable() noexcept {};

  std::vector<BeaconInformationTableEntry>
  load_to_memmory(const dso::datetime<dso::nanoseconds> &t);

}; // BeaconInformationTable

} // namespace dso

#endif