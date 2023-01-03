#include "tides.hpp"
#include <charconv>
#include <cstring>
#include <fstream>
#include <vector>

namespace {
struct OcTideRecordLine {
  int doodson[6];
  int l, m;
  double DelCpl, DelSpl, DelCmi, DelSmi;

  dso::DoodsonOceanTideConstituent as_constituent() const noexcept {
    dso::DoodsonOceanTideConstituent d;
    for (int i=0; i<6; i++) d.doodson[i] = doodson[i];
    d.maxl = d.minl = l;
    d.maxm = m;
    return d;
  }
}; // OcTideRecordLine

bool sameDoodson(const int *i1, const int *i2) noexcept {
  int diff = 0;
  for (int i = 0; i < 6; i++)
    diff += (i1[i] != i2[i]);
  return (diff == 0);
}

int parse_line(const char *line, OcTideRecordLine &rec) noexcept {
  auto skipws = [](const char *str) noexcept -> const char * {
    while (*str && *str == ' ')
      ++str;
    return str;
  };

  // read and resolve Doodson number
  if (dso::doodson2intarray(skipws(line), rec.doodson)) {
    fprintf(stderr,
            "[ERROR] Failed to resolve Doodson number from line: %s "
            "(traceback: %s)\n",
            line, __func__);
    return 1;
  }

  // size of line
  const int sz = std::strlen(line);

  // error status
  int error = 0;

  // Darwin name; start at: 7(=Doodson Number) + 1(whitespace)
  const char *str = line + 7 + 1;
  // Darwin name may have a max size 9in chars) of 4
  // just ignore for now!
  str += 4;

  // read degree and order
  auto cres = std::from_chars(skipws(str), line + sz, rec.l);
  error += (cres.ec != std::errc{});
  str = cres.ptr;
  cres = std::from_chars(skipws(str), line + sz, rec.m);
  error += (cres.ec != std::errc{});
  str = cres.ptr;
  if (error || (rec.l < 0 || rec.m > rec.l)) {
    fprintf(stderr,
            "[ERROR] Failed to resolve SH degree/order from line: %s "
            "(traceback: %s)\n",
            line, __func__);
    return 1;
  }

  // read coefficients in the order: DelC+     DelS+       DelC-     DelS-
  double d[4];
  for (int i = 0; i < 4; i++) {
    cres = std::from_chars(skipws(str), line + sz, d[i]);
    error += (cres.ec != std::errc{});
    str = cres.ptr;
  }
  if (error) {
    fprintf(stderr,
            "[ERROR] Failed to resolve SH coefficients from line: %s "
            "(traceback: %s)\n",
            line, __func__);
    return 1;
  }

  // assign to passed-in instance
  rec.DelCpl = d[0];
  rec.DelSpl = d[1];
  rec.DelCmi = d[2];
  rec.DelSmi = d[3];

  return 0;
}
} // unnamed namespace

constexpr const int MAX_CHARS_IN_OCCOEFFS = 256;

int memmap_octide_coefficients(const char *fn, std::vector<DoodsonOceanTideConstituent> &freqs) noexcept {
}

int dso::inspect_octide_coefficients(
    const char *fn, std::vector<DoodsonOceanTideConstituent> &freqs) noexcept {

  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Failed opening file: %s (traceback: %s)\n", fn,
            __func__);
    return 1;
  }

  freqs.clear();
  char line[MAX_CHARS_IN_OCCOEFFS];

  // read in header comments, first three lines
  for (int i = 0; i < 3; i++)
    fin.getline(line, MAX_CHARS_IN_OCCOEFFS);

  // read in header line, should match!
  const char *HEADER =
      "Doodson Darw  l   m    DelC+     DelS+       DelC-     DelS-";
  const int HSZ = std::strlen(HEADER);
  fin.getline(line, MAX_CHARS_IN_OCCOEFFS);
  if (std::strncmp(HEADER, line, HSZ)) {
    fprintf(stderr, "[ERROR] Unexpected header in file %s! 9traceback: %s)\n",
            fn, __func__);
    fprintf(stderr, "Expected: %s\nFound   : %s\n", HEADER, line);
    return 1;
  }

  // read in coefficients for first constituent
  OcTideRecordLine rec;
  fin.getline(line, MAX_CHARS_IN_OCCOEFFS);
  if (parse_line(line, rec))
    return 1;

  // add this first constituent to the vector
  freqs.emplace_back(rec.as_constituent());

  // keep on reading constituents ...
  int error = 0;
  while (fin.getline(line, MAX_CHARS_IN_OCCOEFFS) && !error) {
    error = parse_line(line, rec);
    if (!error) {
      if (sameDoodson((freqs.end() - 1)->doodson, rec.doodson)) {
        (freqs.end() - 1)->maxl = std::max((freqs.end() - 1)->maxl, rec.l);
        (freqs.end() - 1)->maxm = std::max((freqs.end() - 1)->maxm, rec.m);
      } else {
        freqs.emplace_back(rec.as_constituent());
      }
    }
  }

  // have we read the file through ?
  if (error)
    return error;
  if (!fin.good() && fin.eof())
    return 0;
  return 99;
}
