#ifndef __CONSTEXPR_ORBIT_INTERGATORS_COEFFICIENTS_HPP__
#define __CONSTEXPR_ORBIT_INTERGATORS_COEFFICIENTS_HPP__

#include "gcem.hpp"

namespace dso {

constexpr const double LambdaRKN768 = 1e0 / 20e0;

template <int N, int M> struct RungeKuttaNystromCoefficients {};
template <> struct RungeKuttaNystromCoefficients<4, 3> {
  static constexpr const int S = 4;
  static constexpr const double c[] = {0e0, 1e0 / 4e0, 7e0 / 10e0, 1e0};
  static constexpr const double bhat[] = {1e0 / 14e0, 8e0 / 27e0, 25e0 / 189e0,
                                          0e0};
  static constexpr const double bdothat[] = {1e0 / 14e0, 32e0 / 81e0,
                                             250e0 / 567e0, 5e0 / 54e0};
  static constexpr const double b[] = {-7e0 / 150e0, 67e0 / 150e0, 3e0 / 20e0,
                                       -1e0 / 20e0};
  static constexpr const double bdot[] = {13e0 / 21e0, -20e0 / 27e0,
                                          275e0 / 189e0, -1e0 / 3e0};
  static constexpr const double a[][3] = {
      {0e0, 0e0, 0e0},
      {1e0 / 32e0, 0e0, 0e0},
      {7e0 / 1e3, 119e0 / 500e0, 0e0},
      {1e0 / 14e0, 8e0 / 27e0, 25e0 / 189e0}};
};
template <> struct RungeKuttaNystromCoefficients<6, 4> {
  static constexpr const int S = 6;
  static constexpr const double c[] = {0e0,        1e0 / 10e0,  3e0 / 10e0,
                                       7e0 / 10e0, 17e0 / 25e0, 1e0};
  static constexpr const double bhat[] = {151e0 / 2142e0,    5e0 / 116e0,
                                          385e0 / 1368e0,    55e0 / 168e0,
                                          -6250e0 / 28101e0, 0e0};
  static constexpr const double bdothat[] = {151e0 / 2142e0,      25e0 / 522e0,
                                             275e0 / 684e0,       275e0 / 252e0,
                                             -78125e0 / 112404e0, 1e0 / 12e0};
  static constexpr const double b[] = {
      1349e0 / 157500e0,    7873e0 / 50000e0, 192199e0 / 900000e0,
      521683e0 / 2100000e0, -16e0 / 125e0,    0e0};
  static constexpr const double bdot[] = {
      1349e0 / 157500e0,   7873e0 / 45000e0, 27457e0 / 90000e0,
      521683e0 / 630000e0, -2e0 / 5e0,       1e0 / 12e0};
  static constexpr const double a[][5] = {
      {0e0, 0e0, 0e0, 0e0, 0e0},
      {1e0 / 200e0, 0e0, 0e0, 0e0, 0e0},
      {-1e0 / 2200e0, 1e0 / 22e0, 0e0, 0e0, 0e0},
      {637e0 / 6600e0, -7e0 / 110e0, 7e0 / 33e0, 0e0, 0e0},
      {225437e0 / 1968750e0, -30073e0 / 281250e0, 65569e0 / 281250e0,
       -9367e0 / 984375e0, 0e0},
      {151e0 / 2142e0, 5e0 / 116e0, 385e0 / 1368e0, 55e0 / 168e0,
       -6250e0 / 28101e0}};
};

template <> struct RungeKuttaNystromCoefficients<7, 6> {
  static constexpr const int S = 9;
  static constexpr const double c[] = {0e0,
                                       1e0 / 10e0,
                                       1e0 / 5e0,
                                       3e0 / 8e0,
                                       1e0 / 2e0,
                                       (7e0 - gcem::sqrt(21e0)) / 14e0,
                                       (7e0 + gcem::sqrt(21e0)) / 14e0,
                                       1e0,
                                       1e0};
  static constexpr const double bhat[] = {
      1e0 / 20e0,
      0e0,
      0e0,
      0e0,
      8e0 / 45e0,
      7e0 * (7e0 + gcem::sqrt(21e0)) / 360e0,
      7e0 * (7e0 - gcem::sqrt(21e0)) / 360e0,
      0e0,
      0e0};
  static constexpr const double bdothat[] = {
      1e0 / 20e0,   0e0,          0e0,        0e0, 16e0 / 45e0,
      49e0 / 180e0, 49e0 / 180e0, 1e0 / 20e0, 0e0};
  static constexpr const double b[] = {1e0 / 20e0,
                                       0e0,
                                       0e0,
                                       0e0,
                                       8e0 / 45e0,
                                       7e0 * (7e0 + gcem::sqrt(21e0)) / 360e0,
                                       7e0 * (7e0 - gcem::sqrt(21e0)) / 360e0,
                                       -LambdaRKN768,
                                       LambdaRKN768};
  /*static constexpr const double bdot[] =*/
  static constexpr const double a[][8] = {
      {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
      {1e0 / 200e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
      {1e0 / 150e0, 1e0 / 75e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0},
      {171e0 / 8192e0, 45e0 / 4096e0, 315e0 / 8192e0, 0e0, 0e0, 0e0, 0e0, 0e0},
      {5e0 / 288e0, 25e0 / 528e0, 25e0 / 672e0, 16e0 / 693e0, 0e0, 0e0, 0e0,
       0e0},
      {(1003e0 - 205e0 * gcem::sqrt(21e0)) / 12348e0,
       -25e0 * (751e0 - 173e0 * gcem::sqrt(21e0)) / 90552e0,
       25e0 * (624e0 - 137e0 * gcem::sqrt(21e0)) / 43218e0,
       -128e0 * (361e0 - 79e0 * gcem::sqrt(21e0)) / 237699e0,
       (3411e0 - 745e0 * gcem::sqrt(21e0)) / 24696e0, 0e0, 0e0, 0e0},
      {(793e0 + 187e0 * gcem::sqrt(21e0)) / 12348e0,
       -25e0 * (331e0 + 113e0 * gcem::sqrt(21e0)) / 90552e0,
       25e0 * (1044e0 + 247e0 * gcem::sqrt(21e0)) / 43218e0,
       -128 * (14885e0 + 3779e0 * gcem::sqrt(21e0)) / 9745659e0,
       (3327e0 + 797e0 * gcem::sqrt(21e0)) / 24696e0,
       -(581e0 + 127e0 * gcem::sqrt(21e0)) / 1722e0, 0e0, 0e0},
      {-(157e0 - 3e0 * gcem::sqrt(21e0)) / 378e0,
       25e0 * (143e0 - 10e0 * gcem::sqrt(21e0)) / 2772e0,
       -25e0 * (876e0 + 55e0 * gcem::sqrt(21e0)) / 3969e0,
       1280e0 * (913e0 + 18e0 * gcem::sqrt(21e0)) / 596673e0,
       -(1353e0 + 26e0 * gcem::sqrt(21e0)) / 2268e0,
       7e0 * (1777e0 + 377e0 * gcem::sqrt(21e0)) / 4428e0,
       7e0 * (5e0 - 1e0 * gcem::sqrt(21e0)) / 36e0, 0e0},
      {1e0 / 20e0, 0e0, 0e0, 0e0, 8e0 / 45e0,
       7e0 * (7e0 + 1e0 * gcem::sqrt(21e0)) / 36e0,
       7e0 * (7e0 - 1e0 * gcem::sqrt(21e0)) / 36e0, 0e0}};
};

template <> struct RungeKuttaNystromCoefficients<9,8> {
  static constexpr const int S = 12;
  static constexpr const double c[S] = {
    0e0,
    0.6321437219823755989424902073e-02,
    0.1264287443964751197884980415e-01,
    0.4166666666666666666666666666e+00,
    0.2544604380000000000000000000e-01,
    0.1292344072000000000000000000e+00,
    0.2970774243000000000000000000e+00,
    0.5000000000000000000000000000e+00,
    0.7029225757000000000000000000e+00,
    0.8958905519456625480019182180e+00,
    0.1000000000000000000000000000e+01,
    0.1000000000000000000000000000e+01
  };
  static constexpr const double bhat[S] = {
    -0.3518240655555007629763408588e-02,0e0,0e0,0e0,0.6879703880602338697982077081e-01,0.1177879603832509115294976608e+00,0.1377362395616285338725155827e+00,0.1008501654461988983215135519e+00,0.6135631323226062775059274661e-01,0.1699052322619264917582309574e-01,0e0,0e0
  };
  static constexpr const double bdothat[] = {
    -0.3518240655555007629763408588e-02,
    0e0,
    0e0,
    0e0,
    0.7059336055058270664872682481e-01,
    0.1352694242367759659250257652e+00,
    0.1959479526240353954141973686e+00,
    0.2017003308923977966430271039e+00,
    0.20653307257134660619377379808+00,
    0.1631986677839733565676051585e+00,
    0.3027543199644318023740738961e-01,
    0e0
  };
  static constexpr const double a[][11] = {
    {0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0},
    {0.1998028426208654875176403333e-04, 0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0},
    {0.2664037901611539833568537778e-04, 0.5328075803223079667137075556e-04,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0},
    {0.2865365237083083248442544836e+02, -0.5904090661728807966777976253e+02, 0.3047405980201280273890986972e+02,0e0,0e0,0e0,0e0,0e0,0e0,0e0,0e0},
    {0.1065903941831243410656110659e-03, 0e0, 0.2171588709672456801997590475e-03, 0.1307385389198734629886539400e-08,0e0,0e0,0e0,0e0,0e0,0e0,0e0},
    {0.3030980653858376916419326961e-01, 0e0, -0.7126134845184682180630126302e-01, 0.1568175392509382439415139758e-04, 0.4928662616150566473771384202e-01,0e0,0e0,0e0,0e0,0e0,0e0},
    {-0.1477808541165541247980136312e+00, 0e0, 0.3682723455773522839746690025e+00, 0.4848245662307982600517077975e-03, -0.2155027916119219160929380291e+00, 0.3865397359925407390123095003e-01, 0e0,0e0,0e0,0e0,0e0},
    {-0.1089954028613558927819286465e-01, 0e0, 0e0, 0.7635032467467202877246657217e-02, 0.4547118895489745578797952498e-01, 0.4826544540148205388198081045e-01, 0.3452787346228887673098587201e-01, 0e0,0e0,0e0,0e0},
    {0.5456430942591812422112950282e+00, 0e0, -0.8991472409701229621282353299e+00, -0.1676560070320166624537773864e-01, 0.3402552136626613562636251327e0, 0.1681846222333179210929658770e+00, 0.4235220058882386425502896119e-01, 0.6652778464370135979569806931e-01, 0e0,0e0,0e0},
    {-0.1113406602941552830476082277e+01, 0e0, 0.7546972259090522109972731093e+00, 0.1118600256491523208841810745e-01, 0.7937353844859527762086562061e+00, -0.4656855788249036634796925437e+00, 0.3836619386801684522200183160e+00, -0.2737580302105589161215079526e-01, 0.6449737368017565619020116575e-01,0e0,0e0},
    {0.3114812407323702239935577607e+01,0e0,0e0,0e0,-0.4620858770988909664725067967e+01,0.2485087828585910565589259449e+01,-0.9881138857312002169169839148e+00,0.5265010973962881009401492752e+00,-0.4664170309918856325664075129e-01,0.2921302651339753843370630214e-01,0e0},
    {-0.3518240655555007629763408588e-02,0e0,0e0,0e0,0.6879703880602338697982077081e-01,0.1177879603832509115294976608e+00,0.1377362395616285338725155827e+00,0.1008501654461988983215135519e+00,0.6135631323226062775059274661e-01,0.1699052322619264917582309574e-01,0e0}};
};

/// @struct AdamsBashforthCoefficients
/// @brief Holds Adams-Bashforth Coefficients in the range [0, N+2]
/// Reference: Satellite Orbits: Models, Methods and Applications, Montenbruck
/// et al, 2005; See Eq. 4.56
template <int N> struct AdamsBashforthCoefficients {
  double coeffs[(N > 7) ? (N + 2) : 9];
  static constexpr int num_coeffs() noexcept { return (N > 7) ? (N + 2) : 9; }

  /// @brief Compile-time constructor.
  constexpr AdamsBashforthCoefficients() noexcept : coeffs() {
    constexpr const int num_coeffs = (N > 7) ? (N + 2) : 9;
    coeffs[0] = 1e0;
    coeffs[1] = 1e0 / 2e0;
    coeffs[2] = 5e0 / 12e0;
    coeffs[3] = 3e0 / 8e0;
    coeffs[4] = 251e0 / 720e0;
    coeffs[5] = 95e0 / 288e0;
    coeffs[6] = 19087e0 / 60480e0;
    coeffs[7] = 5257e0 / 17280e0;
    coeffs[8] = 1070017e0 / 3628800e0;
    for (int j = 9; j < num_coeffs; j++) {
      double sum = 0e0;
      for (int k = 0; k < j; k++) {
        sum += coeffs[k] / static_cast<double>(j + 1 - k);
      }
      coeffs[j] = 1e0 - sum;
    }
  }
};

/// @struct AdamsMoultonCoefficients
/// @brief Holds Adams-Moulton Coefficients in the range [0, N+2]
/// Reference: Satellite Orbits: Models, Methods and Applications, Montenbruck
/// et al, 2005; See Eq. 4.64
template <int N> struct AdamsMoultonCoefficients {
  double coeffs[(N > 7) ? (N + 2) : 9];
  static constexpr int num_coeffs() noexcept { return (N > 7) ? (N + 2) : 9; }

  /// @brief Compile-time constructor.
  constexpr AdamsMoultonCoefficients() noexcept : coeffs() {
    constexpr const int num_coeffs = (N > 7) ? (N + 2) : 9;
    coeffs[0] = 1e0;
    coeffs[1] = -1e0 / 2e0;
    coeffs[2] = -1e0 / 12e0;
    coeffs[3] = -1e0 / 24e0;
    coeffs[4] = -19e0 / 720e0;
    coeffs[5] = -3e0 / 160e0;
    coeffs[6] = -863e0 / 60480e0;
    coeffs[7] = -275e0 / 24192e0;
    coeffs[8] = -33953e0 / 3628800e0;
    for (int j = 9; j < num_coeffs; j++) {
      double sum = 0e0;
      for (int k = 0; k < j; k++) {
        sum += coeffs[k] / static_cast<double>(j + 1 - k);
      }
      coeffs[j] = -sum;
    }
  }
};

/// @struct StoermerCowellCoefficients_Delta
/// @brief Holds Stoermer-Cowell Coefficients in the range [0, N+2]
/// Reference: Satellite Orbits: Models, Methods and Applications, Montenbruck
/// et al, 2005; See Eq. 4.90
template <int N> struct StoermerCowellCoefficients_Delta {
  double coeffs[(N > 7) ? (N + 2) : 9];
  static constexpr int num_coeffs() noexcept { return (N > 7) ? (N + 2) : 9; }

  /// @brief Compile-time constructor.
  constexpr StoermerCowellCoefficients_Delta() noexcept : coeffs() {
    constexpr const int num_coeffs = (N > 7) ? (N + 2) : 9;
    coeffs[0] = 1e0;
    coeffs[1] = 0e0;
    coeffs[2] = 1e0 / 12e0;
    coeffs[3] = 1e0 / 12e0;
    coeffs[4] = 19e0 / 240e0;
    coeffs[5] = 3e0 / 40e0;
    coeffs[6] = 863e0 / 12096e0;
    coeffs[7] = 275e0 / 4032e0;
    coeffs[8] = 33953e0 / 518400e0;
    if constexpr (N > 7) {
      constexpr const AdamsMoultonCoefficients<N> gamma;
      for (int j = 9; j < num_coeffs; j++)
        coeffs[j] = (1e0 - j) * gamma.coeffs[j];
    }
  }
};

/// @struct StoermerCowellCoefficients
/// @brief Holds Stoermer-Cowell Coefficients in the range [0, N+2]
/// Reference: Satellite Orbits: Models, Methods and Applications, Montenbruck
/// et al, 2005; See Eq. 4.90
template <int N> struct StoermerCowellCoefficients {
  double coeffs[(N > 7) ? (N + 2) : 9];
  static constexpr int num_coeffs() noexcept { return (N > 7) ? (N + 2) : 9; }

  /// @brief Compile-time constructor.
  constexpr StoermerCowellCoefficients() noexcept : coeffs() {
    constexpr const int num_coeffs = (N > 7) ? (N + 2) : 9;
    coeffs[0] = 1e0;
    coeffs[1] = -1e0;
    coeffs[2] = 1e0 / 12e0;
    coeffs[3] = 0e0;
    coeffs[4] = -1e0 / 240e0;
    coeffs[5] = -1e0 / 240e0;
    coeffs[6] = -221e0 / 60480e0;
    coeffs[7] = -19e0 / 6048;
    coeffs[8] = -9829e0 / 3628800e0;
    if constexpr (N > 7) {
      constexpr const StoermerCowellCoefficients_Delta<N> delta;
      for (int j = 9; j < num_coeffs; j++)
        coeffs[j] = delta.coeffs[j] - delta.coeffs[j - 1];
    }
  }
};

} // namespace dso

#endif
