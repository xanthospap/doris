#include "eop.hpp"
#include "geodesy/units.hpp"
#include "iers2010/iau.hpp"
#include "orbit_integration.hpp"
#include <charconv>
#include <cmath>
#include <cstdio>
#include <datetime/dtfund.hpp>
#include <fstream>
#include <vector>

class RunningStats {
  // see https://www.johndcook.com/blog/standard_deviation/
private:
  int n;
  double m_old, m_new, s_old, s_new;

public:
  void update(double val) noexcept {
    ++n;
    if (n == 1) {
      m_old = m_new = val;
      m_old = 0e0;
    } else {
      m_new = m_old + (val - m_old) / n;
      s_new = s_old + (val - m_old) * (val - m_new);
      m_old = m_new;
      s_old = s_new;
    }
    return;
  }

  double mean() const noexcept { return (n > 0) ? m_new : 0e0; }
  double variance() const noexcept { return ((n > 1) ? m_new / (n - 1) : 0e0); }
  double stddev() const noexcept { return std::sqrt(variance()); }
};

struct Rotary {
  double mjd;
  Eigen::Matrix<double, 3, 3> R;
};

int map_input(const char *fn, std::vector<Rotary> &rots);
inline const char *skipws(const char *line) noexcept {
  const char *str = line;
  while (*str && *str == ' ')
    ++str;
  return str;
}

double gps2tt(double gpst) noexcept {
  constexpr const double offset = (32.184e0 + 19e0) / 86400e0;
  return gpst + offset;
}
double gps2tai(double gpst) noexcept {
  constexpr const double offset = (19e0) / 86400e0;
  return gpst + offset;
}

int main(int argc, char *argv[]) {

  // check input
  if (argc != 3) {
    fprintf(stderr,
            "USAGE: %s [EOP INPUT] [01earthRotation_rotaryMatrix.txt]\n",
            argv[0]);
    return 1;
  }

  // parse reference results
  std::vector<Rotary> refrots;
  if (map_input(argv[2], refrots))
    return 1;

  // fist date in file as datetime instance
  dso::datetime<dso::nanoseconds> d1(
      dso::modified_julian_day(static_cast<int>(refrots[0].mjd)),
      dso::nanoseconds(0));

  // Parse the input EOP data file to create an EopLookUpTable eop_lut
  dso::EopLookUpTable eop_lut;
  const int ref_mjd = d1.as_mjd();
  const int start = ref_mjd - 5;
  const int end = ref_mjd + 6;
  // parse C04 EOPs and convert time-stamps to TT (not UTC)
  if (parse_iers_C04(argv[1], start, end, eop_lut)) {
    fprintf(stderr, "ERROR. Failed collecting EOP data\n");
    return 1;
  }

  // regularize ERP (DUT and LOD)
  eop_lut.regularize();

  // ok, now for every MJD in the reference file
  for (const Rotary &rot : refrots) {

    // IAU 2006/2000A, CIO based, using X,Y series
    Eigen::Matrix<double, 3, 3> rc2i, rpom;
    double era, xlod;
    if (dso::gcrs2itrs(gps2tai(rot.mjd), eop_lut, rc2i, era, rpom, xlod)) {
      fprintf(stderr, "ERROR. Failed to compute matrix!\n");
      return 1;
    }

    // actual matrix:
    const auto rc2ti = Eigen::AngleAxisd(era, -Eigen::Vector3d::UnitZ()) * rc2i;
    const Eigen::Matrix<double, 3, 3> R = (rpom * rc2ti);

    // differences [-]
    printf("%12.5f %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e %+.15e "
           "%+.15e\n",
           rot.mjd, R(0, 0) - rot.R(0, 0), R(0, 1) - rot.R(0, 1),
           R(0, 2) - rot.R(0, 2), R(1, 0) - rot.R(1, 0), R(1, 1) - rot.R(1, 1),
           R(1, 2) - rot.R(1, 2), R(2, 0) - rot.R(2, 0), R(2, 1) - rot.R(2, 1),
           R(2, 2) - rot.R(2, 2));
  }

  return 0;
}

// Read and parse file 01earthRotation_rotartMatrix.txt
// Assume columns:
// gps time [mjd], xx, xy, xz, yx, yy, yz, zx, zy, zz
int map_input(const char *fn, std::vector<Rotary> &rots) {
  constexpr const int MAX_CHARS = 1024;
  char line[MAX_CHARS];
  rots.clear();

  std::ifstream fin(fn);
  if (!fin.is_open()) {
    fprintf(stderr, "ERROR. Failed opening file %s\n", fn);
    return 1;
  }

  // first five lines are GROOPS related
  for (int i=0; i<5;i++) fin.getline(line, MAX_CHARS);

  // read data
  double _data[10];
  while (fin.getline(line, MAX_CHARS)) {
    const char *c = line;
    const int sz = std::strlen(line);
    for (int i = 0; i < 10; i++) {
      auto cres = std::from_chars(skipws(c), line + sz, _data[i]);
      if (cres.ec != std::errc{}) {
        fprintf(stderr, "ERROR Failed resolving line %s\n", line);
        if (!std::strlen(line)) {
          fprintf(stderr, "Line is actually empty, so pretending this never "
                          "happened ...\n");
        } else {
          return 2;
        }
      }
      c = cres.ptr;
    }
    // remember, COLUMN-WISE order!
    rots.push_back(
        {_data[0], Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(_data+1)});
  }

  printf("Number of records collected to compare: %d\n", (int)rots.size());

  if (!fin.good() && fin.eof())
    return 0;
  return 3;
}
