#ifndef __DSO_ODE_FUNCTION_PROTO_HPP__
#define __DSO_ODE_FUNCTION_PROTO_HPP__

#include "eigen3/Eigen/Eigen"
#include "orbit_integration.hpp"

namespace dso {

// Function prototype for first order differential equations
// void f (double x, const Vector y, Vector yp[])
typedef void (*ODEfun)(double x,                 // Independent variable
                       const Eigen::VectorXd &y, // State (function values)
                       /*const Eigen::MatrixXd &Phi,*/
                       Eigen::Ref<Eigen::VectorXd> yp, // Partials/Derivative
                       /*Eigen::MatrixXd &Phip, */
                       dso::IntegrationParameters *params) noexcept;
}// dso

#endif
