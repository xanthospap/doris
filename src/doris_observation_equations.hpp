#ifndef __DSO_DORIS_OBSERVATION_EQUATION_HPP__
#define __DSO_DORIS_OBSERVATION_EQUATION_HPP__

#include "doris_system_info.hpp"

namespace dso {
/// @brief Compute the ionospheric correction of a DORIS 2GHz phase
///        observable in cycles
/// @param[in] L1cycles Observation on L1 carrier (2GHz) [cycles]
/// @param[in] L2cycles Observation on L2 carrier (400MHz) [cycles]
/// @return Correction to be applied on the L1 observable due to ionospheric
///         dealy, in [cycles]. Aka, the corrected (carrier) observable on L1
///         should be L1_iono-free-2GHz = L1cycles + Diono [cycles]
/// @see Lemoine etal, 2016, "Precise orbit determination and station
///         position estimation using DORIS RINEX data", section 2.5.7
inline double carrier_iono_correction(double L1cycles,
                                      double L2cycles) noexcept {
  return (L1cycles - dso::GAMMA_FACTOR_SQRT * L2cycles) / (GAMMA_FACTOR - 1e0);
}

double relativistic_clock_correction(const Eigen::Matrix<double, 3, 1> &recef,
                                     const Eigen::Matrix<double, 3, 1> &vecef,
                                     double GM) noexcept;
double relativistic_clock_correction(const Eigen::Matrix<double, 3, 1> &recef,
                                     const Eigen::Matrix<double, 3, 1> &vecef,
                                     double GM, double J2, double Re) noexcept;

} // namespace dso

#endif
