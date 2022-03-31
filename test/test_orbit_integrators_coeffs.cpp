#include "orbit_integrators_coefficients.hpp"
#include <cassert>
#include <cstdio>
#include <cmath>

constexpr double ABRefCoeffs[] = {
  1e0, 1e0/2e0, 5e0/12e0, 3e0/8e0, 251e0/720e0, 95e0/288e0, 19087e0/60480e0, 5257e0/17280e0, 1070017e0/3628800e0};
constexpr double AMRefCoeffs[] = {
  1e0, -1e0/2e0, -1e0/12e0, -1e0/24e0, -19e0/720e0, -3e0/160e0, -863e0/60480e0, -275e0/24192e0, -33953e0/3628800e0};
constexpr double SDRefCoeffs[] = {
  1e0, 0e0, 1e0/12e0, 1e0/12e0, 19e0/240e0, 3e0/40e0, 863e0/12096e0, 275e0/4032e0, 33953e0/518400e0};
constexpr double SCRefCoeffs[] = {
  1e0, -1e0, 1e0/12e0, 0e0, -1e0/240e0, -1e0/240e0, -221e0/60480e0, -19e0/6048e0, -9829e0/3628800e0};

int main() {
  // test Adams Bashforth Coefficients
  constexpr const dso::AdamsBashforthCoefficients<12> c_ab;
  static_assert(std::abs(c_ab.coeffs[0]-ABRefCoeffs[0])<1e-15);
  static_assert(std::abs(c_ab.coeffs[1]-ABRefCoeffs[1])<1e-15);
  static_assert(std::abs(c_ab.coeffs[2]-ABRefCoeffs[2])<1e-15);
  static_assert(std::abs(c_ab.coeffs[3]-ABRefCoeffs[3])<1e-15);
  static_assert(std::abs(c_ab.coeffs[4]-ABRefCoeffs[4])<1e-15);
  static_assert(std::abs(c_ab.coeffs[5]-ABRefCoeffs[5])<1e-15);
  static_assert(std::abs(c_ab.coeffs[6]-ABRefCoeffs[6])<1e-15);
  static_assert(std::abs(c_ab.coeffs[7]-ABRefCoeffs[7])<1e-15);
  static_assert(std::abs(c_ab.coeffs[8]-ABRefCoeffs[8])<1e-15);
  
  // test Adams Moulton Coefficients
  constexpr const dso::AdamsMoultonCoefficients<12> c_am;
  static_assert(std::abs(c_am.coeffs[0]-AMRefCoeffs[0])<1e-15);
  static_assert(std::abs(c_am.coeffs[1]-AMRefCoeffs[1])<1e-15);
  static_assert(std::abs(c_am.coeffs[2]-AMRefCoeffs[2])<1e-15);
  static_assert(std::abs(c_am.coeffs[3]-AMRefCoeffs[3])<1e-15);
  static_assert(std::abs(c_am.coeffs[4]-AMRefCoeffs[4])<1e-15);
  static_assert(std::abs(c_am.coeffs[5]-AMRefCoeffs[5])<1e-15);
  static_assert(std::abs(c_am.coeffs[6]-AMRefCoeffs[6])<1e-15);
  static_assert(std::abs(c_am.coeffs[7]-AMRefCoeffs[7])<1e-15);
  static_assert(std::abs(c_am.coeffs[8]-AMRefCoeffs[8])<1e-15);
  
  // test StoermerCowellCoefficients_Delta Coefficients
  constexpr const dso::StoermerCowellCoefficients_Delta<12> c_sd;
  static_assert(std::abs(c_sd.coeffs[0]-SDRefCoeffs[0])<1e-15);
  static_assert(std::abs(c_sd.coeffs[1]-SDRefCoeffs[1])<1e-15);
  static_assert(std::abs(c_sd.coeffs[2]-SDRefCoeffs[2])<1e-15);
  static_assert(std::abs(c_sd.coeffs[3]-SDRefCoeffs[3])<1e-15);
  static_assert(std::abs(c_sd.coeffs[4]-SDRefCoeffs[4])<1e-15);
  static_assert(std::abs(c_sd.coeffs[5]-SDRefCoeffs[5])<1e-15);
  static_assert(std::abs(c_sd.coeffs[6]-SDRefCoeffs[6])<1e-15);
  static_assert(std::abs(c_sd.coeffs[7]-SDRefCoeffs[7])<1e-15);
  static_assert(std::abs(c_sd.coeffs[8]-SDRefCoeffs[8])<1e-15);
  
  // test StoermerCowellCoefficients Coefficients
  constexpr const dso::StoermerCowellCoefficients<12> c_sc;
  static_assert(std::abs(c_sc.coeffs[0]-SCRefCoeffs[0])<1e-15);
  static_assert(std::abs(c_sc.coeffs[1]-SCRefCoeffs[1])<1e-15);
  static_assert(std::abs(c_sc.coeffs[2]-SCRefCoeffs[2])<1e-15);
  static_assert(std::abs(c_sc.coeffs[3]-SCRefCoeffs[3])<1e-15);
  static_assert(std::abs(c_sc.coeffs[4]-SCRefCoeffs[4])<1e-15);
  static_assert(std::abs(c_sc.coeffs[5]-SCRefCoeffs[5])<1e-15);
  static_assert(std::abs(c_sc.coeffs[6]-SCRefCoeffs[6])<1e-15);
  static_assert(std::abs(c_sc.coeffs[7]-SCRefCoeffs[7])<1e-15);
  static_assert(std::abs(c_sc.coeffs[8]-SCRefCoeffs[8])<1e-15);

  int sz = sizeof(ABRefCoeffs) / sizeof(double);
  for (int i=0; i<sz; i++)
    printf("Adams-Bashforth Coeff[%2d] %.15f %.15f %.15e\n", i, ABRefCoeffs[i], c_ab.coeffs[i], std::abs(ABRefCoeffs[i] - c_ab.coeffs[i]));
  
  for (int i=sz; i<c_ab.num_coeffs(); i++)
    printf("Adams-Bashforth Coeff[%2d] %17s %.15f %17s\n", i, "--", c_ab.coeffs[i], "--");

  sz = sizeof(AMRefCoeffs) / sizeof(double);
  for (int i=0; i<sz; i++)
    printf("Adams-Moulton Coeff[%2d]   %+.15f %+.15f %.15e\n", i, AMRefCoeffs[i], c_am.coeffs[i], std::abs(AMRefCoeffs[i] - c_am.coeffs[i]));
  
  for (int i=sz; i<c_am.num_coeffs(); i++)
    printf("Adams-Moulton Coeff[%2d]   %18s %+.15f %17s\n", i, "--", c_am.coeffs[i], "--");
  
  sz = sizeof(SDRefCoeffs) / sizeof(double);
  for (int i=0; i<sz; i++)
    printf("Stoermer-Cowell (delta) Coeff[%2d]   %+.15f %+.15f %.15e\n", i, SDRefCoeffs[i], c_sd.coeffs[i], std::abs(SDRefCoeffs[i] - c_sd.coeffs[i]));
  
  for (int i=sz; i<c_sd.num_coeffs(); i++)
    printf("Stoermer-Cowell (delta) Coeff[%2d]   %18s %+.15f %17s\n", i, "--", c_sd.coeffs[i], "--");
  
  sz = sizeof(SCRefCoeffs) / sizeof(double);
  for (int i=0; i<sz; i++)
    printf("Stoermer-Cowell (delta*) Coeff[%2d]   %+.15f %+.15f %.15e\n", i, SCRefCoeffs[i], c_sc.coeffs[i], std::abs(SCRefCoeffs[i] - c_sc.coeffs[i]));
  
  for (int i=sz; i<c_sc.num_coeffs(); i++)
    printf("Stoermer-Cowell (delta*) Coeff[%2d]   %18s %+.15f %17s\n", i, "--", c_sc.coeffs[i], "--");

return 0;
}
