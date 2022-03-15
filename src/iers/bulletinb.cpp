#include "iers_bulletin.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

constexpr const std::size_t MAX_LINE_CHARS = 512;

// retrun a pointer to the first non-whitespace char
inline char *nws(char *line) noexcept {
    char *c = line;
    while (*c && *c == ' ') ++c;
    return c;
}

dso::IersBulletinB::IersBulletinB(const char *fn) {
  assert(std::strlen(fn) < 256);
  std::strcpy(filename, fn);
  stream.open(filename);
  if (!stream.is_open()) {
    fprintf(stderr,
            "[ERROR] Failed opening Bulltin B file %s (traceback: %s)\n",
            filename, __func__);
    throw std::runtime_error("Failed opening file");
  }
}

dso::IersBulletinB::IersBulletinB(dso::IersBulletinB &&other) noexcept {
  std::strcpy(filename, other.filename);
  stream.swap(other.stream);
}

dso::IersBulletinB &
dso::IersBulletinB::operator=(dso::IersBulletinB &&other) noexcept {
  std::strcpy(filename, other.filename);
  stream.swap(other.stream);
  return *this;
}

dso::IersBulletinB::~IersBulletinB() noexcept {
  filename[0] = '\0';
  if (stream.is_open())
    stream.close();
}

int reach_section1(std::ifstream &fin) noexcept {
  char line[MAX_LINE_CHARS], *c, *end;
  fin.seekg(0);

  // first line should be a whitespace chars, and 'BULLETIN B XXX'
  fin.getline(line, MAX_LINE_CHARS);
  c = nws(line);
  if (std::strncmp(c, "BULLETIN B ", 11))
    return 1;
  int bulnum = std::strtol(c+11, &end, 10);
  if (!bulnum || (c == end))
    return 1;

  // ignore next couple of lines ...
  // 1. next line is date [ignored]
  // 2. next is empty [ignored]
  // 3. next is 'Contents are described in ...' [ignored]
  // 4. next is empty [ignored]
  // for (int i = 0; i < 4; i++)
  //   fin.getline(line, MAX_LINE_CHARS);
  // !!BUT!! seems that the above might change and more empty lines may be 
  // injected. Hence, just skip all lines untill we reach
  // '1 - DAILY FINAL VALUES OF x, y, UT1-UTC, dX, dY'
int dummy_it = 0, max_dummy_it = 10;
fin.getline(line, MAX_LINE_CHARS);
while ((std::strncmp(nws(line),
                     "1 - DAILY FINAL VALUES OF x, y, UT1-UTC, dX, dY", 48) &&
        std::strncmp(
            nws(line),
            "1 - DAILY SMOOTHED VALUES OF  x, y, UT1-UTC, UT1R-UTC, dX, dY",
            62)) &&
       (dummy_it++ < max_dummy_it)) {
  fin.getline(line, MAX_LINE_CHARS);
}
if (dummy_it >= max_dummy_it)
  return 500;

// 1 - DAILY FINAL VALUES OF x, y, UT1-UTC, dX, dY (already read, see above)
// Next is:
// 'Angular unit is milliarcsecond (mas), time unit is millisecond (ms).'
fin.getline(line, MAX_LINE_CHARS);
c = nws(line);
if (std::strncmp(c,
                 "Angular unit is milliarcsecond (mas), time unit is "
                 "millisecond ",
                 63))
  return 7;

// 'Upgraded solution from March 1 2017 - consistent with ITRF 2014.'
// well, actually this is non-standard, could be  e.g.
// 'The reference systems are described in the 2006 IERS Annual Report.'
fin.getline(line, MAX_LINE_CHARS);
//c = nws(line);
//if (std::strncmp(
//        c, "Upgraded solution from March 1 2017 - consistent with ITRF 2014.",
//        65))
//  return 8;

// an empty line follows; ignore
fin.getline(line, MAX_LINE_CHARS);

// next two lines describe records; they are standard
fin.getline(line, MAX_LINE_CHARS);
c = nws(line);
if (std::strncmp(c,
                 "DATE     MJD       x       y      UT1-UTC      dX     dY   "
                 "  x err    y err   UT1 err  dX err  dY err",
                 102)){
printf("Expected [%s]\n", "DATE     MJD       x       y      UT1-UTC      dX     dY     x err    y err   UT1 err  dX err  dY err");
printf("Found    [%s]\n", c);
  return 10;
                 }
fin.getline(line, MAX_LINE_CHARS);
c = nws(line);
if (std::strncmp(c,
                 "(0 h UTC)            mas     mas       ms         mas    "
                 "mas     mas      mas      ms     mas     mas",
                 102))
  return 11;

// an empty line follows; ignore
fin.getline(line, MAX_LINE_CHARS);

// next two line should be 'Final values'
fin.getline(line, MAX_LINE_CHARS);
c = nws(line);
if (std::strncmp(c, "Final values", 12))
  return 13;

// skip next three lines ...
for (int i = 0; i < 3; i++)
  fin.getline(line, MAX_LINE_CHARS);
return fin.good();
}

int parse_section1_line(const char *line,
                        dso::IersBulletinB_Section1Block &block) noexcept {
  char *end;
  const char *c = line + 15;
  block.mjd = std::strtol(c, &end, 10);
  if (!block.mjd || c == end)
    return 1;

  double tmp[10];
  for (int i = 0; i < 10; i++) {
    c = end + 1;
    tmp[i] = std::strtod(c, &end);
    if (c == end)
      return i;
  }

  block.x = tmp[0];
  block.y = tmp[1];
  block.dut1 = tmp[2];
  block.dX = tmp[3];
  block.dY = tmp[4];
  block.xerr = tmp[5];
  block.yerr = tmp[6];
  block.dut1err = tmp[7];
  block.dXerr = tmp[8];
  block.dYerr = tmp[9];

  return 0;
}

int dso::IersBulletinB::parse_section1(dso::IersBulletinB_Section1Block *block,
                                       bool include_preliminary) noexcept {
  if (int err_line = reach_section1(stream)) {
    fprintf(stderr,
            "[ERROR] Failed parsing Bulletin B file %s; line nr %d (traceback: "
            "%s)\n",
            filename, err_line, __func__);
    return -1;
  }

  char line[MAX_LINE_CHARS];
  int block_nr = 0;
  // next line to read from stream is a Section1 line (final)
  stream.getline(line, MAX_LINE_CHARS);
  while (*line && line[0] != ' ') {
    if (parse_section1_line(line, block[block_nr])) {
      fprintf(stderr,
              "[ERROR] Failed parsing Bulletin B file %s at Section 1 "
              "(traceback: %s)\n",
              filename, __func__);
      return -2;
    }
    block[block_nr++].type = 'F';
    stream.getline(line, MAX_LINE_CHARS);
  }

  // exit if we do not need the preliminary values
  if (!include_preliminary)
    return block_nr;

  // we should now have read all final values, and reached an empty line
  stream.getline(line, MAX_LINE_CHARS);
  if (std::strncmp(line, " Preliminary extension", 22)) {
    fprintf(stderr,
            "[ERROR] Failed parsing Bulletin B file %s at Section 1 "
            "(traceback: %s)\n",
            filename, __func__);
    return -3;
  }

  // next line to read from stream is a Section1 line (preliminery)
  stream.getline(line, MAX_LINE_CHARS);
  while (*line && line[0] != ' ') {
    if (parse_section1_line(line, block[block_nr])) {
      fprintf(stderr,
              "[ERROR] Failed parsing Bulletin B file %s at Section 1 "
              "(traceback: %s)\n",
              filename, __func__);
      return -4;
    }
    block[block_nr++].type = 'P';
    stream.getline(line, MAX_LINE_CHARS);
  }

  return block_nr;
}

int dso::IersBulletinB::get_section1_at(int imjd,
                                        dso::IersBulletinB_Section1Block &block,
                                        bool include_preliminary) noexcept {
  if (int err_line = reach_section1(stream)) {
    fprintf(stderr,
            "[ERROR] Failed parsing Bulletin B file %s; line nr %d (traceback: "
            "%s)\n",
            filename, err_line, __func__);
    return 1;
  }

  // int imjd = static_cast<int>(std::floor(t.as_mjd()));
  char line[MAX_LINE_CHARS];
  // next line to read from stream is a Section1 line (final)
  stream.getline(line, MAX_LINE_CHARS);
  while (*line && line[0] != ' ') {
    if (parse_section1_line(line, block)) {
      fprintf(stderr,
              "[ERROR] Failed parsing Bulletin B file %s at Section 1 "
              "(traceback: %s)\n",
              filename, __func__);
      return 2;
    }
    block.type = 'F';
    if (block.mjd == imjd)
      return 0;
    stream.getline(line, MAX_LINE_CHARS);
  }

  // exit if we do not need the preliminary values
  if (!include_preliminary)
    return -1;

  // we should now have read all final values, and reached an empty line
  stream.getline(line, MAX_LINE_CHARS);
  if (std::strncmp(line, " Preliminary extension", 22)) {
    fprintf(stderr,
            "[ERROR] Failed parsing Bulletin B file %s at Section 1 "
            "(traceback: %s)\n",
            filename, __func__);
    return 3;
  }

  // next line to read from stream is a Section1 line (preliminery)
  stream.getline(line, MAX_LINE_CHARS);
  while (*line && line[0] != ' ') {
    if (parse_section1_line(line, block)) {
      fprintf(stderr,
              "[ERROR] Failed parsing Bulletin B file %s at Section 1 "
              "(traceback: %s)\n",
              filename, __func__);
      return 4;
    }
    block.type = 'P';
    if (block.mjd == imjd)
      return 0;
    stream.getline(line, MAX_LINE_CHARS);
  }

  return -1;
}