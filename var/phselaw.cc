#include "ggeodesy/units.hpp"
#include <cstdio>

/// Alcatel
struct AlcatelPhaseCorrections_D1 {
  double dnorth{0e0}, deast{0e0}, dup{510e0};
  double plarray[19] = {0.00e0,  2.05e0,  7.24e0,  9.21e0,  6.71e0,  
                      8.14e0, 11.87e0, 12.48e0, 12.28e0, 13.67e0, 
                      13.91e0, 13.01e0, 13.01e0, 11.87e0,  9.70e0, 
                      7.94e0,  4.99e0,  0.41e0, -3.93}; /* 0-90 step=5*/
  constexpr double phase_center_offset_north() const noexcept {
    return dnorth;
  }
  constexpr double phase_center_offset_east() const noexcept {
    return deast;
  }
  constexpr double phase_center_offset_up() const noexcept {
    return dup;
  }
  double get_phase_law_correction(double zenith) const noexcept {
    const double angle(ngpt::rad2deg<double>(zenith));
    if (zenith<0e0 || zenith > ngpt::DPI/2e0) return 0e0;
    int index = std::floor(angle/5e0);
#ifdef DEBUG
    assert(index>=0 && index < sizeof(plarray)-1);
#endif
    return plarray[index] + (angle-5e0*index)*( (plarray[index+1]-plarray[index])/5e0 );
  }
};
/*
struct AlcatelPhaseCorrections_D2 {
  constexpr double dnorth{0e0}, deast{0e0}, dup{335e0};
  constexpr double phase_center_offset_north() const noexcept {
    return dnorth;
  }
  constexpr double phase_center_offset_east() const noexcept {
    return deast;
  }
  constexpr double phase_center_offset_up() const noexcept {
    return dup;
  }
  constexpr double get_phase_law_correction(double zenith) const noexcept {
    return 0e0;
  }
};
struct StarecBCPhaseCorrections_D1 {
  constexpr double dnorth{0e0}, deast{0e0}, dup{487e0};
  constexpr double plarray[] = {0.00e0,  0.06e0, -0.32e0, -1.12e0, -2.87e0, 
                     -4.02e0, -3.44e0, -2.15e0, -1.73e0, -1.73e0, 
                     -0.08e0,  1.37e0,  2.20e0,  5.37e0,  7.02e0, 
                     10.70e0, 13.86e0, 17.27e0, 22.37};
  constexpr double phase_center_offset_north() const noexcept {
    return dnorth;
  }
  constexpr double phase_center_offset_east() const noexcept {
    return deast;
  }
  constexpr double phase_center_offset_up() const noexcept {
    return dup;
  }
};
struct StarecBCPhaseCorrections_D2 {
  constexpr double dnorth{0e0}, deast{0e0}, dup{0e0};
  constexpr double phase_center_offset_north() const noexcept {
    return dnorth;
  }
  constexpr double phase_center_offset_east() const noexcept {
    return deast;
  }
  constexpr double phase_center_offset_up() const noexcept {
    return dup;
  }
  constexpr double get_phase_law_correction(double zenith) const noexcept {
    return 0e0;
  }
};
*/
int main() {
  AlcatelPhaseCorrections_D1 ap;
  for (int z=0; z<91; z++) {
    printf("%10.3f %+15.8f\n", (double)z, ap.get_phase_law_correction(ngpt::deg2rad<double>(z)));
  }

  return 0;
}
