#include <cmath>
#include <algorithm>

struct Jb2008InParams {
    double f10, f10b, s10, s10b, xm10, xm10b, y10, y10b, dstdtc;
};

struct Jb2008OutParams {
    double temp[2];
    double rho;
};

// The alpha are the thermal diffusion coefficients in Eq. (6)
const double alpha[] = {0e0,0e0,0e0,0e0,-0.38e0};

// AL10 is DLOG(10.0)
constexpr const double al10 = 2.3025851e0;

// The AMW are the molecular weights in order: N2, O2, O, Ar, He & H
const double anw[] = {28.0134e0,31.9988e0,15.9994e0,39.9480e0,4.0026e0,1.00797e0};

// AVOGAD is Avogadro's number in mks units (molecules/kmol)
constexpr const double avogad = 02257e26;

// π related consts
constexpr const double PI = iers2010::DPI;
constexpr const double TWOPI = 2e0*iers2010::DPI;
constexpr const double FOURPI = 4e0*iers2010::DPI;
constexpr const double PIOV2 = 2e0/iers2010::DPI;
constexpr const double PIOV4 = 4e0/iers2010::DPI;
constexpr const double DEGRAD = PI / 180e0;

// The FRAC are the assumed sea-level volume fractions in order:
// N2, O2, Ar, and He
const double frac[] = {0.78110e0,0.20955e0,9.3400e-3,1.2890e-5};

// RSTAR is the universal gas-constant in mks units (joules/K/kmol)
constexpr const double rstar = 8314.32e0;

// The R# are values used to establish height step sizes in
// the regimes 90km to 105km, 105km to 500km and 500km upward.
constexpr const double r1 =0.010e0;
constexpr const double r2 =0.025e0;
constexpr const double r3 =0.075e0;

// The WT are weights for the Newton-Cotes Five-Point Quad. formula
const double wt[] = {0.311111111111111e0, 1.422222222222222e0,
                     0.533333333333333e0, 1.422222222222222e0,
                     0.311111111111111e0};

// The CHT are coefficients for high altitude density correction
const double cht[] = {0.22e0,-0.20e-02,0.115e-02,-0.211e-05};


int jb2008(double amjd, const double *sun, const double *sat,
           const Jb2008InParams &in, Jb2008OutParams &out) noexcept {

    double al10n[6], aln[6], tc[4];

    // Equation 14
    double fn = std::pow(in.f10b / 240.e0, 1e0 / 4e0);
    if (fn > 1e0)
      fn = 1e0;
    const double fsb = in.f10b * fn + in.s10b * (1e0 - fn);
    const double tsubc =
        392.4e0 + 3.227e0 * fsb + 0.298e0 * (in.f10 - in.f10b) +
        2.259e0 * (in.s10 - in.s10b) + 0.312e0 * (in.xm10 - in.xm10b) +
        0.178e0 * (in.y10 - in.y10b);
    
    // Equation 15
    const double eta =   0.5e0 * std::abs(sat[1] - sun[1]);
    const double theta = 0.5e0 * std::abs(sat[1] + sun[1]);

    // Equation 16
    const double h = sat[0] - sun[0];
    const double tau = h - 0.64577182e0 + 0.10471976e0 * std::sin(h + 0.75049158e0);
    const double glat  = sat[1];
    const double zht   = sat[2];
    const double glst  = h + PI;
    double glsthr = (glst/DEGRAD)*(24e0/360e0);
    if (glsthr >= 24e0) glsthr = glsthr - 24e0;
    if (glsthr < 0e0) glsthr = glsthr + 24e0;

    // Equation 17
    const double c = std::pow(std::cos(eta),2.5);
    const double s = std::pow(std::sin(theta),2.5);
    const double df = s + (c - s) * std::pow(std::abs(std::cos(0.5e0 * tau)),3e0);
    const double tsubl = tsubc * (1e0 + 0.31e0 * df);

    //  Compute correction to dTc for local solar time and lat correction
    dtsub(in.f10,glsthr,glat,zht,dtclst);

    // Compute the local exospheric temperature.
    // Add geomagnetic storm effect from input dTc value
    out.temp[0] = tsubl + in.dstdtc;
    const double tinf = tsubl + in.dstdtc + dtclst;

    // Equation 9
    const double tsubx = 444.3807e0 + 0.02385e0 * tinf -
                         392.8292e0 * std::exp(-0.0021357e0 * tinf);
    
    // Equation 11
    const double gsubx = 0.054285714e0 * (tsubx - 183e0);

    // The TC array will be an argument in the call to
    // XLOCAL, which evaluates Equation (10) or Equation (13)
    tc[0] = tsubx;
    tc[1] = gsubx;

    // A AND GSUBX/A OF Equation (13)
    tc[2] = (tinf - tsubx)/PIOV2;
    tc[3] = gsubx/tc[2];

    // Equation (5)
    const double z1 = 90e0;
    const double z2 = std::min(sat[2], 105e0);
    const double al = std::log(z2 / z1);
    const double n = idint(al / r1) + 1;
    const double zr = std::exp(al / static_cast<double>(n));
    const double ambar1 = xambar(z1);
    const double tloc1 = xlocal(z1, tc);
    double zend = z1;
    double sum2 = 0e0;
    double ain = ambar1 * xgrav(z1) / tloc1;

    double ambar2,tloc2;
    for (int i = 0; i < n; i++) {
      double z = zend;
      zend = zr * z;
      const double dz = 0.25e0 * (zend - z);
      double sum1 = wt[0] * ain;
      for (int j = 1; j < 5; j++) {
        z = z + dz;
        ambar2 = xambar(z);
        tloc2 = xlocal(z, tc);
        const double gravl = xgrav(z);
        ain = ambar2 * gravl / tloc2;
        sum1 += wt[j] * ain;
      }
      sum2 += dz / sum1;
    }

    const double fact1 = 1000e0/rstar;
    const double rho = 3.46e-6 * ambar2 * tloc1 * std::exp(-fact1*sum2) /ambar1 / tloc2;

    // Equation (2)
    const double anm = avogad * rho;
    const double an  = anm/ambar2;

    // Equation (3)
    const double fact2 = anm/28.960e0;
    aln[0] = std::log(frac[0]*fact2);
    aln[3] = std::log(frac[2]*fact2);
    aln[4] = std::log(frac[3]*fact2);

}