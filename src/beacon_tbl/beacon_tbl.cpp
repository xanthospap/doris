#include "beacon_tbl.hpp"
#include <stdexcept>

std::vector<dso::BeaconInformationTableEntry>
dso::BeaconInformationTable::load_to_memmory(
    const dso::datetime<dso::nanoseconds> &t) {
  std::ifstream fin(fn_);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR. Failed to open file %s (traceback: %s)\n", fn_,
            __func__);
    throw std::runtime_error("Open file error");
  }

  std::vector<BeaconInformationTableEntry> vec;
  vec.reserve(60);

  char line[256], name[4];
  fin.getline(line, 256); // header

  char *end;
  while (fin.getline(line, 256)) {
    std::memcpy(name, line, 4 * sizeof(char));
    const auto from = dso::strptime_ymd_hms<dso::nanoseconds>(line + 5, &end);
    const auto to = dso::strptime_ymd_hms<dso::nanoseconds>(line + 24, &end);
    double height = std::strtod(line + 45, &end);
    if (errno || end == line + 45) {
      fprintf(stderr,
              "ERROR. Failed resolving beacon height from line \'%s\' "
              "(traceback: %s)\n",
              line, __func__);
      throw std::runtime_error("Parse file error");
    }
    if (t >= from && t < to) {
      vec.emplace_back(name, height, from, to);
    }
  }

  return vec;
}