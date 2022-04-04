#ifndef __CONSTEXPR_ORBIT_INTERGATORS_COEFFICIENTS_HPP__
#define __CONSTEXPR_ORBIT_INTERGATORS_COEFFICIENTS_HPP__

namespace dso {

template<int N, int M>
struct RungeKuttaNystromCoefficients {};
template<>
struct RungeKuttaNystromCoefficients<4,3> {
  static constexpr const double c[] = {0e0, 1e0/4e0, 7e0/10e0, 1e0};
  static constexpr const double bhat[] = {1e0/14e0, 8e0/27e0, 25e0/189e0, 0e0};
  static constexpr const double bdothat[] = {1e0/14e0, 32e0/81e0, 250e0/567e0, 5e0/54e0};
  static constexpr const double b[] = {-7e0/150e0, 67e0/150e0, 3e0/20e0, -1e0/20e0};
  static constexpr const double bdot[] = {13e0/21e0, -20e0/27e0, 275e0/189e0, -1e0/3e0};
  static constexpr const double a[][3] = {{0e0,0e0,0e0},{1e0/32e0,0e0,0e0},{7e0/1e3,119e0/500e0,0e0},{1e0/4e0,8e0/27e0,25e0/189e0}};
};
template<>
struct RungeKuttaNystromCoefficients<6,4> {
  static constexpr const double c[] = {0e0, 1e0/10e0, 3e0/10e0, 7e0/10e0, 17e0/25e0, 1e0};
  static constexpr const double bhat[] = {151e0/2142e0, 5e0/116e0, 385e0/1368e0, 55e0/168e0, -6250e0/28101e0,0e0};
  static constexpr const double bdothat[] = {151e0/2142e0, 25e0/522e0, 275e0/684e0, 275e0/252e0, -78125e0/112404e0, 1e0/12e0};
  static constexpr const double b[] = {1349e0 / 157500e0, 7873e0 / 50000e0,
                                       192199e0 / 900000e0, 521683e0 / 2100000e0, -16e0/125e0, 0e0};
  static constexpr const double bdot[] = {1349e0/157500e0, 7873e0/45000e0, 27457e0/90000e0, 521683e0/630000e0, -2e0/5e0, 1e0/12e0};
  static constexpr const double a[][5] = {{0e0,0e0,0e0,0e0,0e0},{1e0/200e0,0e0,0e0,0e0,0e0},{-1e0/2200e0, 1e0/22e0,0e0,0e0,0e0},{637e0/6600e0,-7e0/110e0,7e0/33e0,0e0,0e0},{225437e0/1968750e0,-30073e0/281250e0,65569e0/281250e0,-9367e0/984375e0,0e0}, {151e0/2142e0, 5e0/116e0, 385e0/1368e0,55e0/168e0,-6250e0/28101e0}};
};

/// @struct AdamsBashforthCoefficients
/// @brief Holds Adams-Bashforth Coefficients in the range [0, N+2]
/// Reference: Satellite Orbits: Models, Methods and Applications, Montenbruck
/// et al, 2005; See Eq. 4.56
template<int N>
struct AdamsBashforthCoefficients {
  double coeffs[(N>7)?(N+2):9];
  static constexpr int num_coeffs() noexcept { return (N>7)?(N+2):9; }

  /// @brief Compile-time constructor. 
  constexpr AdamsBashforthCoefficients() noexcept : coeffs() {
    constexpr const int num_coeffs = (N>7)?(N+2):9;
    coeffs[0] = 1e0; 
    coeffs[1] = 1e0/2e0; 
    coeffs[2] = 5e0/12e0; 
    coeffs[3] = 3e0/8e0; 
    coeffs[4] = 251e0/720e0; 
    coeffs[5] = 95e0/288e0; 
    coeffs[6] = 19087e0/60480e0; 
    coeffs[7] = 5257e0/17280e0; 
    coeffs[8] = 1070017e0/3628800e0;
    for (int j=9; j<num_coeffs; j++) {
      double sum = 0e0;
      for (int k=0; k<j; k++) {
        sum += coeffs[k] / static_cast<double>(j+1-k);
      }
      coeffs[j] = 1e0 - sum;
    }
  }
};

/// @struct AdamsMoultonCoefficients
/// @brief Holds Adams-Moulton Coefficients in the range [0, N+2]
/// Reference: Satellite Orbits: Models, Methods and Applications, Montenbruck
/// et al, 2005; See Eq. 4.64
template<int N>
struct AdamsMoultonCoefficients {
  double coeffs[(N>7)?(N+2):9];
  static constexpr int num_coeffs() noexcept { return (N>7)?(N+2):9; }

  /// @brief Compile-time constructor. 
  constexpr AdamsMoultonCoefficients() noexcept : coeffs() {
    constexpr const int num_coeffs = (N>7)?(N+2):9;
    coeffs[0] = 1e0; 
    coeffs[1] = -1e0/2e0; 
    coeffs[2] = -1e0/12e0; 
    coeffs[3] = -1e0/24e0;
    coeffs[4] = -19e0/720e0;
    coeffs[5] = -3e0/160e0;
    coeffs[6] = -863e0/60480e0;
    coeffs[7] = -275e0/24192e0;
    coeffs[8] = -33953e0/3628800e0;
    for (int j=9; j<num_coeffs; j++) {
      double sum = 0e0;
      for (int k=0; k<j; k++) {
        sum += coeffs[k] / static_cast<double>(j+1-k);
      }
      coeffs[j] = -sum;
    }
  }
};

/// @struct StoermerCowellCoefficients_Delta
/// @brief Holds Stoermer-Cowell Coefficients in the range [0, N+2]
/// Reference: Satellite Orbits: Models, Methods and Applications, Montenbruck
/// et al, 2005; See Eq. 4.90
template<int N>
struct StoermerCowellCoefficients_Delta {
  double coeffs[(N>7)?(N+2):9];
  static constexpr int num_coeffs() noexcept { return (N>7)?(N+2):9; }

  /// @brief Compile-time constructor. 
  constexpr StoermerCowellCoefficients_Delta() noexcept : coeffs() {
    constexpr const int num_coeffs = (N>7)?(N+2):9;
    coeffs[0] = 1e0; 
    coeffs[1] = 0e0; 
    coeffs[2] = 1e0/12e0; 
    coeffs[3] = 1e0/12e0;
    coeffs[4] = 19e0/240e0;
    coeffs[5] = 3e0/40e0;
    coeffs[6] = 863e0/12096e0;
    coeffs[7] = 275e0/4032e0;
    coeffs[8] = 33953e0/518400e0;
    if constexpr (N>7) {
      constexpr const AdamsMoultonCoefficients<N> gamma;
      for (int j=9; j<num_coeffs; j++)
        coeffs[j] = (1e0 - j) * gamma.coeffs[j];
    }
  }
};

/// @struct StoermerCowellCoefficients
/// @brief Holds Stoermer-Cowell Coefficients in the range [0, N+2]
/// Reference: Satellite Orbits: Models, Methods and Applications, Montenbruck
/// et al, 2005; See Eq. 4.90
template<int N>
struct StoermerCowellCoefficients {
  double coeffs[(N>7)?(N+2):9];
  static constexpr int num_coeffs() noexcept { return (N>7)?(N+2):9; }

  /// @brief Compile-time constructor. 
  constexpr StoermerCowellCoefficients() noexcept : coeffs() {
    constexpr const int num_coeffs = (N>7)?(N+2):9;
    coeffs[0] = 1e0; 
    coeffs[1] = -1e0;
    coeffs[2] = 1e0/12e0;
    coeffs[3] = 0e0;
    coeffs[4] = -1e0/240e0;
    coeffs[5] = -1e0/240e0;
    coeffs[6] = -221e0/60480e0;
    coeffs[7] = -19e0/6048;
    coeffs[8] = -9829e0/3628800e0;
    if constexpr (N>7) {
      constexpr const StoermerCowellCoefficients_Delta<N> delta;
      for (int j=9; j<num_coeffs; j++)
        coeffs[j] = delta.coeffs[j] - delta.coeffs[j-1];
    }
  }
};

}// dso

#endif
