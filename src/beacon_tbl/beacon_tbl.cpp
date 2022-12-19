#include "beacon_tbl.hpp"
#include <charconv>
#include <stdexcept>

namespace {
inline const char *skipws(const char *line) noexcept {
  const char *c = line;
  while (*c && *c == ' ')
    ++c;
  return c;
}
} // unnamed namespace

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

  while (fin.getline(line, 256)) {
    std::memcpy(name, line, 4 * sizeof(char));
    const auto from =
        dso::strptime_ymd_hms<dso::nanoseconds>(skipws(line + 5), nullptr);
    const auto to =
        dso::strptime_ymd_hms<dso::nanoseconds>(skipws(line + 24), nullptr);
    double height = 0e0;
    auto cres = std::from_chars(skipws(line + 45), line + 256, height);
    if (cres.ec != std::errc{}) {
      fprintf(stderr,
              "[ERROR]. Failed resolving beacon height from line \'%s\' "
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
