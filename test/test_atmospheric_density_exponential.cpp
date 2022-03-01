#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "atmosphere.hpp"
#include "iers2010/iersc.hpp"
#include "doctest.h"

using dso::air_density_models::exponential::density;

/// Test from Vallado, Chapter 8.6.2, Example 8-4

TEST_CASE("Testing the Exponential Model for Atmospheric Density.") {
    double sat_radius = 7125.3489e3; // [m]
    double sat_altitude = (sat_radius - iers2010::Re)*1e-3; //[km]
    double rho = density(sat_altitude);
    REQUIRE(std::abs(rho-2.1219854e-14) < 1e-15);
}