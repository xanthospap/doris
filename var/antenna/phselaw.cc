#include "ggeodesy/units.hpp"
#include <cstdio>
#include <iostream>

enum class GroundAntennaType : int { Alcatel, Starec_B, Starec_C };

template<GroundAntennaType T, int Freq> struct AntennaOffsetTraits { };

template<>
struct AntennaOffsetTraits<GroundAntennaType::Alcatel, 1> {
  static constexpr double offset[3] = { 0e0, 0e0, 510e0 };
  static constexpr double plarray[] = {0.00e0,  2.05e0,  7.24e0,  9.21e0,  6.71e0,  
                      8.14e0, 11.87e0, 12.48e0, 12.28e0, 13.67e0, 
                      13.91e0, 13.01e0, 13.01e0, 11.87e0,  9.70e0, 
                      7.94e0,  4.99e0,  0.41e0, -3.93}; /* 0-90 step=5*/
  static constexpr double pcv(double zenith, int& out_of_bounds) noexcept {
    const double angle = ngpt::rad2deg<double>(zenith);
    out_of_bounds = 0;
    if (zenith<0e0 || zenith > ngpt::DPI/2e0) {
      return 0e0;
      out_of_bounds=1;
    }
    int index = std::floor(angle/5e0);
#ifdef DEBUG
    assert(index>=0 && index < sizeof(plarray)-1);
#endif
    return plarray[index] + (angle-5e0*index)*( (plarray[index+1]-plarray[index])/5e0 );
  }
};

template<>
struct AntennaOffsetTraits<GroundAntennaType::Alcatel, 2> {
  static constexpr double offset[3] = { 0e0, 0e0, 335e0 };
  static constexpr double pcv(double zenith, int& out_of_bounds) noexcept {
    out_of_bounds = 0;
    return 0e0;
  }
};

template<>
struct AntennaOffsetTraits<GroundAntennaType::Starec_B, 1> {
  static constexpr double offset[3] = {0e0, 0e0, 487e0};
  static constexpr double plarray[] = {0.00e0,  0.06e0, -0.32e0, -1.12e0, -2.87e0, 
                     -4.02e0, -3.44e0, -2.15e0, -1.73e0, -1.73e0, 
                     -0.08e0,  1.37e0,  2.20e0,  5.37e0,  7.02e0, 
                    10.70e0, 13.86e0, 17.27e0, 22.37};/* 0-90 step=5*/
  static constexpr double pcv(double zenith, int& out_of_bounds) noexcept {
    const double angle(ngpt::rad2deg<double>(zenith));
    out_of_bounds = 0;
    if (zenith<0e0 || zenith > ngpt::DPI/2e0) {
      out_of_bounds = 1;
      return 0e0;
    }
    int index = std::floor(angle/5e0);
#ifdef DEBUG
    assert(index>=0 && index < sizeof(plarray)-1);
#endif
    return plarray[index] + (angle-5e0*index)*( (plarray[index+1]-plarray[index])/5e0 );
  }
};

template<>
struct AntennaOffsetTraits<GroundAntennaType::Starec_B, 2> {
  static double constexpr offset[3] = {0e0, 0e0, 0e0};
  static constexpr double pcv(double zenith, int& out_of_bounds) noexcept {
    out_of_bounds = 0;
    return 0e0;
  }
};

template<>
struct AntennaOffsetTraits<GroundAntennaType::Starec_C, 1> {
  static constexpr double offset[3] = {0e0, 0e0, 487e0};
  static constexpr double plarray[] = {0.00e0,  0.06e0, -0.32e0, -1.12e0, -2.87e0, 
                     -4.02e0, -3.44e0, -2.15e0, -1.73e0, -1.73e0, 
                     -0.08e0,  1.37e0,  2.20e0,  5.37e0,  7.02e0, 
                    10.70e0, 13.86e0, 17.27e0, 22.37};/* 0-90 step=5*/
  static constexpr double pcv(double zenith, int& out_of_bounds) noexcept {
    const double angle(ngpt::rad2deg<double>(zenith));
    out_of_bounds = 0;
    if (zenith<0e0 || zenith > ngpt::DPI/2e0) {
      out_of_bounds = 1;
      return 0e0;
    }
    int index = std::floor(angle/5e0);
#ifdef DEBUG
    assert(index>=0 && index < sizeof(plarray)-1);
#endif
    return plarray[index] + (angle-5e0*index)*( (plarray[index+1]-plarray[index])/5e0 );
  }
};

template<>
struct AntennaOffsetTraits<GroundAntennaType::Starec_C, 2> {
  static constexpr double offset[3] = {0e0, 0e0, 0e0};
  static constexpr double pcv(double zenith, int& out_of_bounds) noexcept {
    out_of_bounds = 0;
    return 0e0;
  }
};

template<GroundAntennaType T, int Freq> 
struct AntennaOffset : AntennaOffsetTraits<T, Freq> {
  static constexpr double dnorth() noexcept { 
    return AntennaOffsetTraits<T, Freq>::offset[0];
  }
  static constexpr double deast() noexcept { 
    return AntennaOffsetTraits<T, Freq>::offset[1];
  }
  static constexpr double dup() noexcept { 
    return AntennaOffsetTraits<T, Freq>::offset[2];
  }
  static constexpr double pcv(double zenithdeg, int& out_of_bounds) noexcept {
    return AntennaOffsetTraits<T, Freq>::pcv(zenithdeg, out_of_bounds);
  }
};

int main() {
  static_assert(AntennaOffset<GroundAntennaType::Starec_B, 1>::dnorth() == 0e0);
  // AntennaOffset<GroundAntennaType::Starec_B, 5>::dnorth();
  static_assert(AntennaOffset<GroundAntennaType::Starec_B, 1>::dnorth() == 
    AntennaOffset<GroundAntennaType::Starec_C, 1>::dnorth());

  int d1_out_of_bounds, d2_out_of_bounds;
  printf("Zenith      Alcatel_1  Alcatel_2  Starec_1  Starec_2\n");
  for (double zenith=0e0; zenith <= 90e0; zenith += 0.3e0) {
    double zrad = ngpt::deg2rad<double>(zenith);
    double alc1 = AntennaOffset<GroundAntennaType::Alcatel, 1>::pcv(zrad, d1_out_of_bounds);
    double alc2 = AntennaOffset<GroundAntennaType::Alcatel, 2>::pcv(zrad, d2_out_of_bounds);
    double str1 = AntennaOffset<GroundAntennaType::Starec_B, 1>::pcv(zrad, d1_out_of_bounds);
    double str2 = AntennaOffset<GroundAntennaType::Starec_B, 2>::pcv(zrad, d2_out_of_bounds);
    printf("%+10.4f %+10.4f %+10.4f %+10.4f %+10.4f\n", zenith, alc1, alc2, str1, str2);
    if (d1_out_of_bounds || d2_out_of_bounds) {
      std::cerr << "[ERROR] Out-of-bounds zenith angle!";
      return 2;
    }
  }
  return 0;
}
