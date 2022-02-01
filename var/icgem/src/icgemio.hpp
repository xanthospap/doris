#ifndef __ICGEM_POTENTIAL_IO_HPP__
#define __ICGEM_POTENTIAL_IO_HPP__

#include <fstream>
#include <string>
#include "harmonic_coeffs.hpp"

class Icgem {
public:
  typedef std::ifstream::pos_type pos_type;

private:
    std::string filename;
    pos_type data_section_pos{0};
    char product_type[64] = {'\0'};
    char modelname[128] = {'\0'};
    char tide_system[64] = "unknown";
    char norm[64] = "fully_normalized";
    char errors[64] = {'\0'};
    double earth_gravity_constant{0e0};
    double radius{0e0};
    int max_degree{0};

public:
    Icgem(const char *fn) noexcept : filename(fn) {};
    Icgem(const Icgem&) noexcept = delete;
    Icgem& operator=(const Icgem&) noexcept = delete;
    ~Icgem() noexcept {};

    int degree() const noexcept { return max_degree; }

    #ifdef DEBUG
    void print_details();
    #endif

    int parse_header() noexcept;
    int parse_data(int l, int m, HarmonicCoeffs *coeffs) noexcept;
}; //Icgem

#endif