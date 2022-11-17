
///< Nominal values of solid Earth tide external potential Love numbers.
///< IERS2010, Table 6.3
struct {
    int n,m;
    double knm, knm_Rplus, Reknm, Imknm, knm_Iplus;
} const LoveK[] = {{2, 0, 0.29525e0, -0.00087e0, 0.30190e0, 0e0, -0.00089e0},
             {2, 1, 0.29470e0, -0.00079e0, 0.29830e0, -0.00144e0, -0.00080e0},
             {2, 2, 0.29801e0, -0.00057e0, 0.30102e0, -0.00130e0, -0.00057e0},
             {3, 0, 0.093e0, 0e0, 0e0, 0e0, 0e0},
             {3,1,0.093e0, 0e0, 0e0, 0e0, 0e0},
             {3,2,0.093e0, 0e0, 0e0, 0e0, 0e0},
             {3,3,0.094e0, 0e0, 0e0, 0e0, 0e0}};

struct SolidEarthTide {
    const double GM, GM_moon, GM_sun;
    const double Re;
    dso::AssociatedLegendreFunctions p;

    int solid_earth_tide_step1(double Rmoon, double Rsun, double mlon, double slon) noexcept {
      const double RRm = Re/Rmoon;
      const double RRm3 = RRm *RRm *RRm;
      const double RRs = Re / Rsun;
      const double RRs3 = RRs * RRs * RRs;
      const double GMme = GM_moon / GM;
      const double GMse = GM_sun / GM;
      const double sml = std::sin(mlon);
      const double ssl = std::sin(slon);
      const double cml = std::cos(mlon);
      const double csl = std::cos(slon);

      // n = 2, m = 0
      double fac = (LoveK[0].knm) / (LoveK[0].n * 2 + 1);
      const double dc20_moon = fac * (GMme) * (RRm3) * p(2, 0);
      const double dc20_sun = fac *  (GMse) * (RRs3) * p(2, 0);

      // n = 2, m = 1
      fac = (LoveK[1].knm) / (LoveK[1].n * 2 + 1);
      const double dc21_moon = fac * (GMme) * (RRm3) * p(2,1) * cml;
      const double ds21_moon = fac * (GMme) * (RRm3) * p(2,1) * sml;
      const double dc21_sun =  fac * (GMse) * (RRs3) * p(2,1) * csl;
      const double ds21_sun =  fac * (GMse) * (RRs3) * p(2,1) * ssl;

      // n = 2, m = 2
      fac = (LoveK[2].knm) / (LoveK[2].n * 2 + 1);
      const double dc22_moon = fac * (GMme) * (RRm3) * p(2,2) * std::cos(2e0*mlon);
      const double ds22_moon = fac * (GMme) * (RRm3) * p(2,2) * std::sin(2e0*mlon);
      const double dc22_sun =  fac * (GMse) * (RRs3) * p(2,2) * std::cos(2e0*slon);
      const double ds22_sun =  fac * (GMse) * (RRs3) * p(2,2) * std::sin(2e0*slon);

      // n = 3, m = 0
      const double RRm4 = RRm3 * RRm;
      const double RRs4 = RRs3 * RRs;
      fac = (LoveK[3].knm) / (LoveK[3].n * 2 + 1);
      const double dc30_moon = fac * (GMme) * (RRm4) * p(3,0);
      const double dc30_sun  = fac * (GMse) * (RRs4) * p(3,0);
      
      // n = 3, m = 1
      fac = (LoveK[4].knm)/(LoveK[4].n*2+1);
      const double dc31_moon = fac * (GMme) * (RRm4) * p(3,1) * cml;
      const double ds31_moon = fac * (GMme) * (RRm4) * p(3,1) * sml;
      const double dc31_sun =  fac * (GMse) * (RRs4) * p(3,1) * csl;
      const double ds31_sun =  fac * (GMse) * (RRs4) * p(3,1) * ssl;

      // n = 3, m = 2
      fac = (LoveK[5].knm) / (LoveK[5].n * 2 + 1);
      const double dc32_moon = fac * (GMme) * (RRm4) * p(3,2) * std::cos(2e0*mlon);
      const double ds32_moon = fac * (GMme) * (RRm4) * p(3,2) * std::sin(2e0*mlon);
      const double dc32_sun =  fac * (GMse) * (RRs4) * p(3,2) * std::cos(2e0*slon);
      const double ds32_sun =  fac * (GMse) * (RRs4) * p(3,2) * std::sin(2e0*slon);
      
      // n = 3, m = 3
      fac = (LoveK[6].knm) / (LoveK[6].n * 2 + 1);
      const double dc33_moon = fac * (GMme) * (RRm4) * p(3,3) * std::cos(3e0*mlon);
      const double ds33_moon = fac * (GMme) * (RRm4) * p(3,3) * std::sin(3e0*mlon);
      const double dc33_sun =  fac * (GMse) * (RRs4) * p(3,3) * std::cos(3e0*slon);
      const double ds33_sun =  fac * (GMse) * (RRs4) * p(3,3) * std::sin(3e0*slon);
      
      // n = 4, m = 0
      fac = (LoveK[0].knm_Rplus) / 5;
      const double dc40_moon = fac * (GMme) * (RRm3) * p(2,0);
      const double dc40_sun  = fac * (GMse) * (RRs3) * p(2,0);
      
      // n = 4, m = 1
      fac = (LoveK[1].knm_Rplus) / 5e0;
      const double dc41_moon = fac * (GMme) * (RRm3) * p(2,1) * cml;
      const double ds41_moon = fac * (GMme) * (RRm3) * p(2,1) * sml;
      const double dc41_sun =  fac * (GMse) * (RRs3) * p(2,1) * csl;
      const double ds41_sun =  fac * (GMse) * (RRs3) * p(2,1) * ssl;

      // n = 4, m = 2
      fac = (LoveK[2].knm_Rplus) / 5e0;
      const double dc42_moon = fac * (GMme) * (RRm3) * p(2,2) * std::cos(2e0*mlon);
      const double ds42_moon = fac * (GMme) * (RRm3) * p(2,2) * std::sin(2e0*mlon);
      const double dc42_sun =  fac * (GMse) * (RRs3) * p(2,2) * std::cos(2e0*slon);
      const double ds42_sun =  fac * (GMse) * (RRs3) * p(2,2) * std::sin(2e0*slon);
      
      return 0;
    }

int operator()(double mjd_tai) noexcept 
{
    p.compute(Moon_lat);
    const double Rmoon;

    p.compute(Sun_lat);
}
}; // SolidEarthTide
