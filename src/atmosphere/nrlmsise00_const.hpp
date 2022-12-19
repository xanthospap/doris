#ifndef __DSO_CVERS_NRLMSISE00_CONST_HPP__
#define __DSO_CVERS_NRLMSISE00_CONST_HPP__

namespace dso::nrlmsise00::detail {

// Constants (taken from physics_constants.f90)

// Boltzman Constant - CODATA 2010
constexpr const double boltzk = 1.380648813e0;

// Avogadro Constant - International Avogadro cohordination 2011
constexpr const double navgdr = 6.0221407818e0;

// Gas constant in J/K/mol
constexpr const double rgasmol = boltzk * navgdr;

// Molecular weight of water
constexpr const double watmolwgt = 18.01528e0; // g/mol

// Mean dry air molecular weight
constexpr const double airmolwgt = 28.96443e0; // g/mol

// Ratio of mean molecular weight of water to that of dry air
constexpr const double wgtfrac = watmolwgt / airmolwgt;

// Gas constant for dry air in J/K/kg
constexpr const double rgas = (rgasmol / airmolwgt) * 1e3;

// Gas constant for water in J/K/kg
constexpr const double rwat = (rgasmol / watmolwgt) * 1e3;

// 0 C in Kelvin
constexpr const double tzero = 273.16e0;

// Standard Gravity (m/sec**2) 3rd CGPM
constexpr const double egrav = 9.80665e0;

// Earth radius in meters
constexpr const double earthrad = 6.371229e+06;
constexpr const double erkm = earthrad / 1e3;

// Angular velocity of rotation of Earth
constexpr const double eomeg = 7.2921159e-05;
constexpr const double eomeg2 = 2e0 * eomeg;

// Hydrostatic coefficient
constexpr const double gmr = egrav * airmolwgt / rgasmol;

// Specific heat at constant pressure for dry air J/kg/K
constexpr const double cpd = 1005.46e0;

// Specific heat at constant pressure for moist air J/kg/K
constexpr const double cpv = 1869.46e0;

// Specific heat of water at 15 Celsius J/kg/K
constexpr const double cpw = 4186.95e0;

// Specific heat of water at 0 Celsius J/kg/K
constexpr const double cpw0 = 4218.0e0;

// Derived
constexpr const double rgovrw = rgas / rwat;
constexpr const double rwovrg = rwat / rgas;
constexpr const double rgovcp = rgas / cpd;
constexpr const double rgovg = rgas / egrav;
constexpr const double govrg = egrav / rgas;

constexpr const double regrav = 1e0 / egrav;
constexpr const double rrgas = 1e0 / rgas;
constexpr const double rcpd = 1e0 / cpd;

constexpr const double r100gas = rgasmol * 100e0;
constexpr const double nearzero = 0.000001e0;
constexpr const double dr = 1.72142e-2;
constexpr const int mn1 = 5;
constexpr const int mn2 = 4;
constexpr const int mn3 = 5;

// const double zn1[] = {120e0, 110e0, 100e0, 90e0, 72.e0};
// const double zn2[] = {72.5e0, 55e0, 45e0, 32.e0};
// const double zn3[] = {32.5e0, 20e0, 15e0, 10e0, 0e0};
// const double mt[] = {48, 0, 4, 16, 28, 32, 40, 1, 49, 14, 17};
// const double altl[] = {200e0, 300e0, 160e0, 250e0, 240e0, 450e0, 320e0,
// 450e0}; const double alpha[] = {-0.3e0, 0e0, 0e0, 0e0, 0.1e0, 0e0, -0.3e0,
// 0e0, 0e0};

/* PARMB */
// static double gsurf;
// static double re;
//
///* GTS3C */
// static double dd;
//
///* DMIX */
// static double dm04, dm16, dm28, dm32, dm40, dm01, dm14;
//
///* MESO7 */
// static double tn1[5];
// static double tn2[4];
// static double tn3[5];
// static double tgn1[2];
// static double tgn2[2];
// static double tgn3[2];
//
///* LPOLY */
// static double dfa;
// static double plg[4][9];
// static double ctloc, stloc;
// static double c2tloc, s2tloc;
// static double s3tloc, c3tloc;
// static double apdf, apt[4];

} // namespace dso::nrlmsise00::detail

#endif
