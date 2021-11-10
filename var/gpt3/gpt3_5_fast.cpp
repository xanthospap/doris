#include <cmath>
#include "ggdatetime/dtcalendar.hpp"
#include "ggeodesy/units.hpp"

// function [p,T,dT,Tm,e,ah,aw,la,undu,Gn_h,Ge_h,Gn_w,Ge_w] = gpt3_5_fast (mjd,lat,lon,h_ell,it,grid)
int gpt3_5_fast(const dso::datetime<dso::nanoseconds> &t, const double *lat, const double *lon, const double *hell, int it, const char *grid_file, int num_stations) noexcept {
  
  // determine the GPT3 coefficients
  // mean gravity in m/s**2
  constexpr double gm = 9.80665e0;
  // molar mass of dry air in kg/mol
  constexpr double dMtr = 28.965e-3;
  // universal gas constant in J/K/mol
  constexpr double Rg = 8.3143e0;

  // get t's doy
  dso::ydoy_date ydoy (t.as_ydoy());
  int idoy = ydoy.__doy.as_underlying_type();
    double fdoy = static_cast<double>(idoy) + t.fractional_days();

  double cosfy=0e0, coshy=0e0, sinfy=0e0, sinhy=0e0;
  // factors for amplitudes
  if (it != 1) {
      cosfy = std::cos(fdoy/365.25e0*2*pi);  // coefficient for A1
      coshy = std::cos(fdoy/365.25e0*4*pi);  // coefficient for B1
      sinfy = std::sin(fdoy/365.25e0*2*pi);  // coefficient for A2
      sinhy = std::sin(fdoy/365.25e0*4*pi);  // coefficient for B2
  }

  // used for indexing lines
  int index[4];
  int bilinear=0;

  // loop over stations
  for (int k=0; k<num_stations; k++) {
    
    // only positive longitude in degrees
    double plon = dso::rad2deg(lon[k] + (lon[k]<0)*2e0*pi);
    // transform to polar distance in degrees
    double ppod = dso::rad2deg(-lat[k] + pi/2e0);

    // find the index (line in the grid file) of the nearest point
    int ipod = std::floor((ppod+5)/5);
    int ilon = std::floor((plon+5)/5);

    // normalized (to one) differences, can be positive or negative
    int diffpod = (ppod - (ipod*5 - 2.5))/5;
    int difflon = (plon - (ilon*5 - 2.5))/5;
    if (ipod == 37)
        ipod = 36;
    if (ilon == 73)
        ilon = 1;
    else if (ilon == 0)
        ilon = 72;

    // get the number of the corresponding line
    indx[0] = (ipod - 1)*72 + ilon;
    
    // near the poles: nearest neighbour interpolation, otherwise: bilinear
    bilinear = 0;
    if (ppod > 2.5 && ppod < 177.5)
           bilinear = 1;
    
    // case of nearest neighbourhood
    if (bilinear == 0) {

        int ix = indx[0];
        
        // transforming ellipsoidal height to orthometric height
        undu[k] = u_grid[ix];
        double hgt = h_ell[k]-undu[k];

        // pressure, temperature at the height of the grid
        double T0 = T_grid[ix,0] + T_grid[ix,1]*cosfy + T_grid[ix,2]*sinfy + T_grid[ix,3]*coshy + T_grid[ix,4]*sinhy;
        double p0 = p_grid[ix,0] + p_grid[ix,2]*cosfy + p_grid[ix,2]*sinfy + p_grid[ix,3]*coshy + p_grid[ix,4]*sinhy;
         
        // specific humidity
        double Q = Q_grid[ix,0] + Q_grid[ix,1]*cosfy + Q_grid[ix,2]*sinfy + Q_grid[ix,3]*coshy + Q_grid[ix,4]*sinhy;
            
        // lapse rate of the temperature
        dT(k) = dT_grid[ix,0] + dT_grid[ix,1]*cosfy + dT_grid[ix,2]*sinfy + dT_grid[ix,3]*coshy + dT_grid[ix,4]*sinhy; 
        
        // temperature lapse rate in degrees / km
        dT(k) = dT[k]*1000;

        // station height - grid height
        double redh = hgt - Hs_grid[ix];

        // temperature at station height in Celsius
        T(k) = T0 + dT[k]*redh - 273.15e0;

        // virtual temperature in Kelvin
        Tv = T0*(1e0+0.6077e0*Q);
        
        double c = gm*dMtr/(Rg*Tv);
        
        // pressure in hPa
        p(k) = (p0*std::exp(-c*redh))/100e0;
            
        // hydrostatic and wet coefficients ah and aw 
        ah(k) = ah_grid[ix,0] + ah_grid[ix,1]*cosfy + ah_grid[ix,2]*sinfy + ah_grid[ix,3]*coshy + ah_grid[ix,4]*sinhy;
        aw(k) = aw_grid[ix,0] + aw_grid[ix,1]*cosfy + aw_grid[ix,2]*sinfy + aw_grid[ix,3]*coshy + aw_grid[ix,4]*sinhy;
    
    	// water vapour decrease factor la
        la(k) = la_grid[ix,0] +
                la_grid[ix,1]*cosfy + la_grid[ix,2]*sinfy +
                la_grid[ix,3]*coshy + la_grid[ix,4]*sinhy; 
		
		// mean temperature Tm
        Tm(k) = Tm_grid[ix,0] +
                Tm_grid[ix,1]*cosfy + Tm_grid[ix,2]*sinfy +
                Tm_grid[ix,3]*coshy + Tm_grid[ix,4]*sinhy;
        
        // north and east gradients (total, hydrostatic and wet)
        Gn_h(k) = Gn_h_grid[ix,0] + Gn_h_grid[ix,1]*cosfy + Gn_h_grid[ix,2]*sinfy + Gn_h_grid[ix,3]*coshy + Gn_h_grid[ix,4]*sinhy;
        Ge_h(k) = Ge_h_grid[ix,0] + Ge_h_grid[ix,1]*cosfy + Ge_h_grid[ix,2]*sinfy + Ge_h_grid[ix,3]*coshy + Ge_h_grid[ix,4]*sinhy;
        Gn_w(k) = Gn_w_grid[ix,0] + Gn_w_grid[ix,1]*cosfy + Gn_w_grid[ix,2]*sinfy + Gn_w_grid[ix,3]*coshy + Gn_w_grid[ix,4]*sinhy;
        Ge_w(k) = Ge_w_grid[ix,0] + Ge_w_grid[ix,1]*cosfy + Ge_w_grid[ix,2]*sinfy + Ge_w_grid[ix,3]*coshy + Ge_w_grid[ix,4]*sinhy;
		
		// water vapor pressure in hPa
		double e0 = Q*p0/(0.622e0+0.378e0*Q)/100e0; // on the grid
		e(k) = e0*(100e0*p(k)/p0)^(la(k)+1e0); // on the station height - (14) Askne and Nordius, 1987
    
    } else { //% bilinear interpolation
        
        ipod1 = ipod + sign(diffpod);
        ilon1 = ilon + sign(difflon);
        if ilon1 == 73
            ilon1 = 1;
        end
        if ilon1 == 0
            ilon1 = 72;
        end
        
        % get the number of the line
        indx(2) = (ipod1 - 1)*72 + ilon;  % along same longitude
        indx(3) = (ipod  - 1)*72 + ilon1; % along same polar distance
        indx(4) = (ipod1 - 1)*72 + ilon1; % diagonal
}

