#ifndef __DSO_ODE_FUNCTION_PROTO_HPP__
#define __DSO_ODE_FUNCTION_PROTO_HPP__

#include "eigen3/Eigen/Eigen"
#include "orbit_integration.hpp"

namespace dso {

/* Function prototype for first order differential equations */
typedef void (*ODEfun)(double x,                 // Independent variable
                       const Eigen::VectorXd &y, // State (function values)
                       Eigen::Ref<Eigen::VectorXd> yp, // Partials/Derivative
                       dso::IntegrationParameters *params) noexcept;
} /* namespace dso */

#endif
