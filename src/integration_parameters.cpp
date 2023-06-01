#include "orbit_integration.hpp"
#include <stdexcept>
dso::IntegrationParameters::IntegrationParameters(const dso::EopLookUpTable &eoptable_,
                        dso::EarthGravity *earth_gravity,
                        dso::OceanTide *ocean_tide,
                        dso::SolidEarthTide *solid_earth_tide,
                        dso::SolidEarthPoleTide *solid_earth_pole_tide,
                        dso::OceanPoleTide *ocean_pole_tide,
                        const char *pck_kernel,
                        const char *dtm2020datafile)
      : eopLUT(eoptable_), 
        egravity(earth_gravity),
        octide(ocean_tide),
        setide(solid_earth_tide),
        psetide(solid_earth_pole_tide),
        poctide(ocean_pole_tide),
        Dtm20(dtm2020datafile),
        geopotential(egravity->max_degree()),
        V(egravity->max_degree() + 3, egravity->max_degree() + 3),
        W(egravity->max_degree() + 3, egravity->max_degree() + 3) {
    /* Sun and Moon gravitational constants in SI units (m^3/s^2) */
    assert(!dso::get_sun_moon_GM(pck_kernel, GMSun, GMMon, true));
    /* make sure that workspace is ok */
    int max_degree = egravity->max_degree();
    int max_order = egravity->max_order();
    if (octide) {
      if ((octide->max_degree() > egravity->max_degree()) || (octide->max_order() > egravity->max_order())) {
        fprintf(stderr, "[WRNNG] Max degree/order for Ocean Tides is larger than for gravity! (traceback: %s)\n", __func__);
        max_degree = std::max(max_degree, octide->max_degree());
        max_order = std::max(max_order, octide->max_order());
      }
    }
    if (setide) {
      if ((setide->max_degree() > egravity->max_degree()) || (setide->max_order() > egravity->max_order())) {
        fprintf(stderr, "[WRNNG] Max degree/order for Solid Earth Tides is larger than for gravity! (traceback: %s)\n", __func__);
        max_degree = std::max(max_degree, setide->max_degree());
        max_order = std::max(max_order, setide->max_order());
      }
    }
    if (poctide) {
      if ((poctide->max_degree() > egravity->max_degree()) || (poctide->max_order() > egravity->max_order())) {
        fprintf(stderr, "[WRNNG] Max degree/order for Ocean Pole Tides is larger than for gravity! (traceback: %s)\n", __func__);
        max_degree = std::max(max_degree, poctide->max_degree());
        max_order = std::max(max_order, poctide->max_order());
      }
    }
    if (psetide) {
      if ((psetide->max_degree() > egravity->max_degree()) || (psetide->max_order() > egravity->max_order())) {
        fprintf(stderr, "[WRNNG] Max degree/order for Solid Earth Pole Tides is larger than for gravity! (traceback: %s)\n", __func__);
        max_degree = std::max(max_degree, psetide->max_degree());
        max_order = std::max(max_order, psetide->max_order());
      }
    }
  }
