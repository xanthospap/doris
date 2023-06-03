#ifndef __DSO_SAT_OCCULTATION_HPP__
#define __DSO_SAT_OCCULTATION_HPP__

#include "eigen3/Eigen/Eigen"
#include "iers2010/iersc.hpp"
#include <cmath>

namespace dso {

/* @brief Compute shadow factor based on the conic shadow model
 *
 * We only consider here the occultation due to Earth; Moon's shadowing is not 
 * considered at al.
 * The model is described in : Zhang et al, 2019, "Study of satellite shadow
 * function model considering the overlapping parts of Earth shadow and Moon
 * shadow and its application to GPS satellite orbit determination"
 *
 * @param[in] rsat Satellite position in ECI [m]
 * @param[in] rsun Sun position in ECI [m]
 * @param[in] Rs Radius of the Sun [m]
 * @param[in] Re Radius of the Earth [m]
 */
double conic_shadow_factor(const Eigen::Matrix<double, 3, 1> &rsat,
                    const Eigen::Matrix<double, 3, 1> &rsun, double Rs,
                    double Re = iers2010::Re) noexcept;
} /* namespace dso */

#endif
