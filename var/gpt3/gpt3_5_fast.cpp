#include <cmath>
#include "gpt3.hpp"
#include "ggdatetime/dtcalendar.hpp"
#include "ggeodesy/units.hpp"

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

// function [p,T,dT,Tm,e,ah,aw,la,undu,Gn_h,Ge_h,Gn_w,Ge_w] = gpt3_5_fast (mjd,lat,lon,h_ell,it,grid)
int dso::gpt3_5_fast(const dso::datetime<dso::nanoseconds> &t, const double *lat,
                const double *lon, const double *hell, int num_stations, int it,
                const char *grid_file) noexcept {

  // determine the GPT3 coefficients
  // mean gravity in m/s**2
  constexpr double gm = 9.80665e0;
  // molar mass of dry air in kg/mol
  constexpr double dMtr = 28.965e-3;
  // universal gas constant in J/K/mol
  constexpr double Rg = 8.3143e0;

  // get t's doy and fraction of day
  dso::ydoy_date ydoy(t.as_ydoy());
  int idoy = ydoy.__doy.as_underlying_type();
  double fdoy = static_cast<double>(idoy) + t.fractional_days();

  double cosfy = 0e0, coshy = 0e0, sinfy = 0e0, sinhy = 0e0;
  // factors for amplitudes
  if (it != 1) {
    cosfy = std::cos(fdoy / 365.25e0 * 2 * pi); // coefficient for A1
    coshy = std::cos(fdoy / 365.25e0 * 4 * pi); // coefficient for B1
    sinfy = std::sin(fdoy / 365.25e0 * 2 * pi); // coefficient for A2
    sinhy = std::sin(fdoy / 365.25e0 * 4 * pi); // coefficient for B2
  }

  // used for indexing lines
  int indx[4];
  int bilinear = 0;

  // loop over stations
  for (int k = 0; k < num_stations; k++) {

    // only positive longitude in degrees
    double plon = dso::rad2deg(lon[k] + (lon[k] < 0) * 2e0 * pi);
    // transform to polar distance in degrees
    double ppod = dso::rad2deg(-lat[k] + pi / 2e0);

    // find the index (line in the grid file) of the nearest point
    int ipod = std::floor((ppod + 5) / 5);
    int ilon = std::floor((plon + 5) / 5);

    // normalized (to one) differences, can be positive or negative
    int diffpod = (ppod - (ipod * 5 - 2.5)) / 5;
    int difflon = (plon - (ilon * 5 - 2.5)) / 5;
    if (ipod == 37)
      ipod = 36;
    if (ilon == 73)
      ilon = 1;
    else if (ilon == 0)
      ilon = 72;

    // get the number of the corresponding line
    indx[0] = (ipod - 1) * 72 + ilon;

    // near the poles: nearest neighbour interpolation, otherwise: bilinear
    bilinear = 0;
    if (ppod > 2.5 && ppod < 177.5)
      bilinear = 1;

    // case of nearest neighbourhood
    if (bilinear == 0) {

      int ix = indx[0];

      // transforming ellipsoidal height to orthometric height
      undu[k] = u_grid[ix];
      double hgt = h_ell[k] - undu[k];

      // pressure, temperature at the height of the grid
      double T0 = T_grid[ix, 0] + T_grid[ix, 1] * cosfy +
                  T_grid[ix, 2] * sinfy + T_grid[ix, 3] * coshy +
                  T_grid[ix, 4] * sinhy;
      double p0 = p_grid[ix, 0] + p_grid[ix, 2] * cosfy +
                  p_grid[ix, 2] * sinfy + p_grid[ix, 3] * coshy +
                  p_grid[ix, 4] * sinhy;

      // specific humidity
      double Q = Q_grid[ix, 0] + Q_grid[ix, 1] * cosfy + Q_grid[ix, 2] * sinfy +
                 Q_grid[ix, 3] * coshy + Q_grid[ix, 4] * sinhy;

      // lapse rate of the temperature
      dT[k] = dT_grid[ix, 0] + dT_grid[ix, 1] * cosfy + dT_grid[ix, 2] * sinfy +
              dT_grid[ix, 3] * coshy + dT_grid[ix, 4] * sinhy;

      // temperature lapse rate in degrees / km
      dT[k] = dT[k] * 1000e0;

      // station height - grid height
      double redh = hgt - Hs_grid[ix];

      // temperature at station height in Celsius
      T[k] = T0 + dT[k] * redh - 273.15e0;

      // virtual temperature in Kelvin
      double Tv = T0 * (1e0 + 0.6077e0 * Q);

      double c = gm * dMtr / (Rg * Tv);

      // pressure in hPa
      p[k] = (p0 * std::exp(-c * redh)) / 100e0;

      // hydrostatic and wet coefficients ah and aw
      ah[k] = ah_grid[ix, 0] + ah_grid[ix, 1] * cosfy + ah_grid[ix, 2] * sinfy +
              ah_grid[ix, 3] * coshy + ah_grid[ix, 4] * sinhy;
      aw[k] = aw_grid[ix, 0] + aw_grid[ix, 1] * cosfy + aw_grid[ix, 2] * sinfy +
              aw_grid[ix, 3] * coshy + aw_grid[ix, 4] * sinhy;

      // water vapour decrease factor la
      la[k] = la_grid[ix, 0] + la_grid[ix, 1] * cosfy + la_grid[ix, 2] * sinfy +
              la_grid[ix, 3] * coshy + la_grid[ix, 4] * sinhy;

      // mean temperature Tm
      Tm[k] = Tm_grid[ix, 0] + Tm_grid[ix, 1] * cosfy + Tm_grid[ix, 2] * sinfy +
              Tm_grid[ix, 3] * coshy + Tm_grid[ix, 4] * sinhy;

      // north and east gradients (total, hydrostatic and wet)
      Gn_h[k] = Gn_h_grid[ix, 0] + Gn_h_grid[ix, 1] * cosfy +
                Gn_h_grid[ix, 2] * sinfy + Gn_h_grid[ix, 3] * coshy +
                Gn_h_grid[ix, 4] * sinhy;
      Ge_h[k] = Ge_h_grid[ix, 0] + Ge_h_grid[ix, 1] * cosfy +
                Ge_h_grid[ix, 2] * sinfy + Ge_h_grid[ix, 3] * coshy +
                Ge_h_grid[ix, 4] * sinhy;
      Gn_w[k] = Gn_w_grid[ix, 0] + Gn_w_grid[ix, 1] * cosfy +
                Gn_w_grid[ix, 2] * sinfy + Gn_w_grid[ix, 3] * coshy +
                Gn_w_grid[ix, 4] * sinhy;
      Ge_w[k] = Ge_w_grid[ix, 0] + Ge_w_grid[ix, 1] * cosfy +
                Ge_w_grid[ix, 2] * sinfy + Ge_w_grid[ix, 3] * coshy +
                Ge_w_grid[ix, 4] * sinhy;

      // water vapor pressure in hPa
      double e0 = Q * p0 / (0.622e0 + 0.378e0 * Q) / 100e0; // on the grid
      e[k] =
          e0 * (100e0 * p(k) / p0) ^
          (la(k) + 1e0); // on the station height - (14) Askne and Nordius, 1987

    } else { //% bilinear interpolation

      int ipod1 = ipod + sgn(diffpod);
      int ilon1 = ilon + sgn(difflon);
      if (ilon1 == 73)
        ilon1 = 1;
      else if (ilon1 == 0)
        ilon1 = 72;

      // get the number of the line
      indx[1] = (ipod1 - 1) * 72 + ilon;  // along same longitude
      indx[2] = (ipod - 1) * 72 + ilon1;  // along same polar distance
      indx[3] = (ipod1 - 1) * 72 + ilon1; // diagonal

      // transforming ellipsoidal height to orthometric height : 
      // Hortho = -N + Hell
      double undul[4], hgt[4];
      grid.undul_slice(indx, undul);
      for (int i=0; i<4; i++) {
        hgt[i] = h_ell[k] - undul[i];
      }

      // pressure, temperature at the height of the grid
      double T0[4], p0[4];
      for (int i = 0; i < 4; i++) {
        T0[i] = T_grid[indx[i], 0] + T_grid[indx[i], 1] * cosfy +
                T_grid[indx[i], 2] * sinfy + T_grid[indx[i], 3] * coshy +
                T_grid[indx[i], 4] * sinhy;
      }
      for (int i = 0; i < 4; i++) {
        p0[i] = p_grid[indx[i], 0] + p_grid[indx[i], 1] * cosfy +
                p_grid[indx[i], 2] * sinfy + p_grid[indx[i], 3] * coshy +
                p_grid[indx[i], 4] * sinhy;
      }
      // humidity
      double Ql[4];
      for (int i=0; i<4; i++) {
      Ql[i] = Q_grid[indx[i], 0] + Q_grid[indx[i], 1] * cosfy +
                  Q_grid[indx[i], 2] * sinfy + Q_grid[indx[i], 3] * coshy +
                  Q_grid[indx[i], 4] * sinhy;
      }

      // reduction = stationheight - gridheight
      //double Hs1 = Hs_grid[indx];
      //double redh = hgt - Hs1;
      double Hs1[4], redh[4];
      grid.Hsgrid_slice(index, Hs1);
      for (int i=0; i<4; i++) redh[i] = hgt[i] - Hs1[i]

      // lapse rate of the temperature in degree / m
      double dT1[4];
      for (int i=0; i<4; i++) {
      dTl[i] = dT_grid[indx[i], 0] + dT_grid[indx[i], 1] * cosfy +
                   dT_grid[indx[i], 2] * sinfy + dT_grid[indx[i], 3] * coshy +
                   dT_grid[indx[i], 4] * sinhy;
      }

      // temperature reduction to station height
      double Tl[4];
      for (int i=0; i<4; i++) {
        Tl[i] = T0[i] + dTl[i]*redh[i] - 273.15e0;
      }

      // virtual temperature
      double Tv[4], c[4];
      for (int i = 0; i < 4; i++) {
        Tv[i] = T0[i] * (1e0 + 0.6077e0 * Ql[i]);
        c[i] = gm * dMtr / (Rg * Tv[i]);
      }

      // pressure in hPa
      double pl[4];
      for (int i=0; i<4; i++) {
      pl[i] = (p0[i]*std::exp(-c[i]*redh[i])) / 100e0;
      }

      // hydrostatic and wet coefficients ah and aw
      ahl = ah_grid[indx, 0] + ah_grid[indx, 1] * cosfy +
            ah_grid[indx, 2] * sinfy + ah_grid[indx, 3] * coshy +
            ah_grid[indx, 4] * sinhy;
      awl = aw_grid[indx, 0] + aw_grid[indx, 1] * cosfy +
            aw_grid[indx, 2] * sinfy + aw_grid[indx, 3] * coshy +
            aw_grid[indx, 4] * sinhy;

      // water vapour decrease factor la
      lal = la_grid[indx, 0] + la_grid[indx, 1] * cosfy +
            la_grid[indx, 2] * sinfy + la_grid[indx, 3] * coshy +
            la_grid[indx, 4] * sinhy;

      // mean temperature of the water vapor Tm
      Tml = Tm_grid[indx, 0] + Tm_grid[indx, 1] * cosfy +
            Tm_grid[indx, 2] * sinfy + Tm_grid[indx, 3] * coshy +
            Tm_grid[indx, 4] * sinhy;

      // north and east gradients(total, hydrostatic and wet)
      Gn_hl = Gn_h_grid[indx, 0] + Gn_h_grid[indx, 1] * cosfy +
              Gn_h_grid[indx, 2] * sinfy + Gn_h_grid[indx, 3] * coshy +
              Gn_h_grid[indx, 4] * sinhy;
      Ge_hl = Ge_h_grid[indx, 0] + Ge_h_grid[indx, 1] * cosfy +
              Ge_h_grid[indx, 2] * sinfy + Ge_h_grid[indx, 3] * coshy +
              Ge_h_grid[indx, 4] * sinhy;
      Gn_wl = Gn_w_grid[indx, 0] + Gn_w_grid[indx, 1] * cosfy +
              Gn_w_grid[indx, 2] * sinfy + Gn_w_grid[indx, 3] * coshy +
              Gn_w_grid[indx, 4] * sinhy;
      Ge_wl = Ge_w_grid[indx, 0] + Ge_w_grid[indx, 1] * cosfy +
              Ge_w_grid[indx, 2] * sinfy + Ge_w_grid[indx, 3] * coshy +
              Ge_w_grid[indx, 4] * sinhy;

      // water vapor pressure in hPa
      double e0 = Ql.*p0./ (0.622e0 + 0.378e0 * Ql) / 100e0; // on the grid
      double el =
          e0.*(100. * pl./ p0).^
          (lal + 1); // on the station height - (14)Askne and Nordius, 1987

      int dnpod1 = std::abs(diffpod); // distance nearer point
      int dnpod2 = 1 - dnpod1;        // distance to distant point
      int dnlon1 = std::abs(difflon);
      int dnlon2 = 1 - dnlon1;

      // pressure
      double R1 = dnpod2 * pl(1) + dnpod1 * pl(2);
      double R2 = dnpod2 * pl(3) + dnpod1 * pl(4);
      p[k] = dnlon2 * R1 + dnlon1 * R2;

      // temperature
      R1 = dnpod2 * Tl(1) + dnpod1 * Tl(2);
      R2 = dnpod2 * Tl(3) + dnpod1 * Tl(4);
      T[k] = dnlon2 * R1 + dnlon1 * R2;

      // temperature in degree per km
      R1 = dnpod2 * dTl(1) + dnpod1 * dTl(2);
      R2 = dnpod2 * dTl(3) + dnpod1 * dTl(4);
      dT[k] = (dnlon2 * R1 + dnlon1 * R2) * 1000;

      // water vapor pressure in hPa
      R1 = dnpod2 * el(1) + dnpod1 * el(2);
      R2 = dnpod2 * el(3) + dnpod1 * el(4);
      e[k] = dnlon2 * R1 + dnlon1 * R2;

      // ah and aw
      R1 = dnpod2 * ahl(1) + dnpod1 * ahl(2);
      R2 = dnpod2 * ahl(3) + dnpod1 * ahl(4);
      ah[k] = dnlon2 * R1 + dnlon1 * R2;
      R1 = dnpod2 * awl(1) + dnpod1 * awl(2);
      R2 = dnpod2 * awl(3) + dnpod1 * awl(4);
      aw[k] = dnlon2 * R1 + dnlon1 * R2;

      // undulation
      R1 = dnpod2 * undul(1) + dnpod1 * undul(2);
      R2 = dnpod2 * undul(3) + dnpod1 * undul(4);
      undu[k] = dnlon2 * R1 + dnlon1 * R2;

      // water vapor decrease factor la
      R1 = dnpod2 * lal(1) + dnpod1 * lal(2);
      R2 = dnpod2 * lal(3) + dnpod1 * lal(4);
      la[k] = dnlon2 * R1 + dnlon1 * R2;

      // gradients
      R1 = dnpod2 * Gn_hl(1) + dnpod1 * Gn_hl(2);
      R2 = dnpod2 * Gn_hl(3) + dnpod1 * Gn_hl(4);
      Gn_h[k] = (dnlon2 * R1 + dnlon1 * R2);
      R1 = dnpod2 * Ge_hl(1) + dnpod1 * Ge_hl(2);
      R2 = dnpod2 * Ge_hl(3) + dnpod1 * Ge_hl(4);
      Ge_h[k] = (dnlon2 * R1 + dnlon1 * R2);
      R1 = dnpod2 * Gn_wl(1) + dnpod1 * Gn_wl(2);
      R2 = dnpod2 * Gn_wl(3) + dnpod1 * Gn_wl(4);
      Gn_w[k] = (dnlon2 * R1 + dnlon1 * R2);
      R1 = dnpod2 * Ge_wl(1) + dnpod1 * Ge_wl(2);
      R2 = dnpod2 * Ge_wl(3) + dnpod1 * Ge_wl(4);
      Ge_w[k] = (dnlon2 * R1 + dnlon1 * R2);

      // mean temperature of the water vapor Tm
      R1 = dnpod2 * Tml(1) + dnpod1 * Tml(2);
      R2 = dnpod2 * Tml(3) + dnpod1 * Tml(4);
      Tm[k] = dnlon2 * R1 + dnlon1 * R2;
    }
  }

  return 0;
}
