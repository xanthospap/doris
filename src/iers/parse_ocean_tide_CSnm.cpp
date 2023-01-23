#include "tides.hpp"
#include <algorithm>
#include <charconv>
#include <cstring>
#include <exception>
#include <fstream>
#include <iers2010/doodson.hpp>
#include <vector>

constexpr const int MAX_CHARS_IN_OCCOEFFS = 256;

namespace {
struct OcTideRecordLine {
  dso::DoodsonNumber doodson;
  int l, m;
  double DelCpl, DelSpl, DelCmi, DelSmi;

  dso::DoodsonOceanTideConstituent as_constituent() const noexcept {
    return dso::DoodsonOceanTideConstituent(doodson, l, m);
  }
}; // OcTideRecordLine

int parse_doodson_iers_str(const char *str, dso::DoodsonNumber &d) noexcept {
  // check that the input string has size >=7 and that all but the 3rd char
  // are actually numeric values
  for (int i = 0; i < 7; i++) {
    if (!str[i] || !(std::isdigit(*(unsigned char *)(str + i)) || i == 3)) {
      fprintf(stderr,
              "[ERROR] Failed to resolve Doodson number from string %s "
              "(traceback: %s)\n",
              str, __func__);
      return 1;
    }
  }

  // first three ints, for variables τ, s, h
  d(0) = str[0] - '0';
  d(1) = str[1] - '0' - 5;
  d(2) = str[2] - '0' - 5;

  // 4th character should be either a '.' or a ','
  if (str[3] != '.' && str[3] != ',') {
    fprintf(stderr,
            "[ERROR] Failed to resolve Doodson number from string %s "
            "(traceback: %s)\n",
            str, __func__);
    return 1;
  }

  // next three characters, for variables p, N', ps
  d(3) = str[4] - '0' - 5;
  d(4) = str[5] - '0' - 5;
  d(5) = str[6] - '0' - 5;

  // all done
  return 0;
}

int parse_line(const char *line, OcTideRecordLine &rec) noexcept {
  auto skipws = [](const char *str) noexcept -> const char * {
    while (*str && *str == ' ')
      ++str;
    return str;
  };

  // read and resolve Doodson number. IERS-2010 uses +5 for the Doodson number
  if (parse_doodson_iers_str(skipws(line), rec.doodson)) {
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

/// @note The input argument vec cannot be declared const, cause then the
///       return type of the find_if call will be an iterator to const,
///       hence different from guess
auto find_doodson_in_vec(
    std::vector<dso::DoodsonOceanTideConstituent> &vec,
    const dso::DoodsonNumber &d,
    std::vector<dso::DoodsonOceanTideConstituent>::iterator guess) noexcept {

  if (guess != vec.end() && guess->doodson_number() == d) [[likely]]
    return guess;

  auto it = std::find_if(
      vec.begin(), vec.end(),
      [D = std::cref(d)](const dso::DoodsonOceanTideConstituent &element) {
        return (element.doodson_number() == D);
      });
  /*[=](const dso::DoodsonOceanTideConstituent &el) {
    return el.doodson_number() == d;
  });*/

  return it;
}

} // unnamed namespace

int dso::memmap_octide_coefficients(
    const char *fn, std::vector<dso::DoodsonOceanTideConstituent> &freqs,
    int max_degree, int max_order, int num_header_lines, double scale) noexcept {
  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "[ERROR] Failed opening file: %s (traceback: %s)\n", fn,
            __func__);
    return 1;
  }

  freqs.clear();
  char line[MAX_CHARS_IN_OCCOEFFS];

  // read in header comments, first three lines
  for (int i = 0; i < num_header_lines; i++)
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

  // keep on reading constituents ...
  OcTideRecordLine rec;
  auto it = freqs.end();
  int error = 0;
  while (fin.getline(line, MAX_CHARS_IN_OCCOEFFS) && !error) {
    error = parse_line(line, rec);
    if (!error) {
      // find constintuent in vector
      it = find_doodson_in_vec(freqs, rec.doodson, it);
      if (it != freqs.cend()) {
        // constituent found, add elements for given degree and order
        // (if needed)
        if (rec.l <= max_degree && rec.m <= max_order) {
          it->delCp(rec.l, rec.m) = rec.DelCpl*scale;
          it->delSp(rec.l, rec.m) = rec.DelSpl*scale;
          it->delCm(rec.l, rec.m) = rec.DelCmi*scale;
          it->delSm(rec.l, rec.m) = rec.DelSmi*scale;
        }
      } else {
        // else, add new constituent
        freqs.emplace_back(rec.doodson, max_degree, max_order);
        // freqs.push_back(dso::DoodsonOceanTideConstituent(rec.doodson,
        // max_degree, max_order));
        it = std::prev(freqs.end());
        it->clear_coefficients();
        // it->print_matrix_sizes();
        if (rec.l <= max_degree && rec.m <= max_order) {
          it->delCp(rec.l, rec.m) = rec.DelCpl*scale;
          it->delSp(rec.l, rec.m) = rec.DelSpl*scale;
          it->delCm(rec.l, rec.m) = rec.DelCmi*scale;
          it->delSm(rec.l, rec.m) = rec.DelSmi*scale;
        }
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
