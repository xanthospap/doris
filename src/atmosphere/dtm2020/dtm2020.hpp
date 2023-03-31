#ifndef __DSO_DTM2020_DRAG_HPP__
#define __DSO_DTM2020_DRAG_HPP__

#include "datetime/dtcalendar.hpp"
#include "solar_activity_data.hpp"

namespace dso {

namespace details {

double DIN__day_of_year(const dso::TwoPartDate &t) noexcept;

/* TODO is this correct ? */
double DIN__local_hours_radians(const dso::TwoPartDate &t, double lon) noexcept;

/* conversion geographic to geomagnetic latitude&longitude */
int geogm(double xlat, double xlong, double &glat, double &glon) noexcept;

class Dtm2020Input {
//private:
public: /* TODO only for debuging */
  double lon; /* longitude [rad] ITRF */
  double lat; /* latitude [rad] ITRF */
  double alt; /* altitude [km] */
  double lh;  /* local hours [rad] */
  double doy; /* day of yesr [1-366] plus fractional part (?) */

public:
  double altitude_km() const noexcept { return alt; }
  double longitude() const noexcept { return lon; }
  double latitude() const noexcept { return lat; }
  double local_hours_radians() const noexcept { return lh; }
  double day_of_year() const noexcept { return doy; }

  bool validate_params() const noexcept {
    return (lh >= 0e0 && lh <= 24e0 && doy >= 0e0 && doy < 366e0 &&
            alt > 120e0);
  }

  /* !WARNING!
   * 1. assumes lon has already been set
   * 2. does not perform parameter check (use validate_params())
   * For how to properly use this function, see Dtm2020::set()
   */
  void set_date(const dso::TwoPartDate &t) noexcept {
    lh = DIN__local_hours_radians(t, lon);
    doy = DIN__day_of_year(t);
  }

  /* !WARNING!
   * 1. does not perform parameter check (use validate_params())
   * 2. You should now call set_date() to compute local time/hours based on
   *    the new longitude
   * For how to properly use this function, see Dtm2020::set()
   * @param[in] longitude Longitude [-π, π] [rad]
   * @param[in] latitude  Latitude [-π/2, π/2] [rad]
   * @param[in] latitude  Altitude [km]. Should be > 120 km
   */
  void set_position(double longitude, double latitude,
                    double altitude) noexcept {
    lon = longitude;
    lat = latitude;
    alt = altitude;
  }

}; // Dtm2020Input

class Dtm2020Output {
  double d[10];
  /*
   *   0: partial density of atomic hydrogen (in gram/cm3)
   *   1: partial density of helium
   *   2: partial density of atomic oxygen
   *   3: partial density of molecular nitrogen
   *   4: partial density of molecular oxygen
   *   5: partial density of atomic nitrogen
   * ---+-------------------------------------------------
   *   6: [tz] temperature at altitude
   *   7: [tinf] exospheric temperature
   *   8: [ro] total density (in gram/cm3)
   *   9: [wmm] mean molecular mass (in gram)
   */
public:
  double hydrogen() const noexcept { return d[0]; }
  double helium() const noexcept { return d[1]; }
  double atomic_oxygen() const noexcept { return d[2]; }
  double molecular_nitrogen() const noexcept { return d[3]; }
  double molecular_oxygen() const noexcept { return d[4]; }
  double atomic_nitrogen() const noexcept { return d[5]; }
  double temperature() const noexcept { return d[6]; }
  double exospheric_temperature() const noexcept { return d[7]; }
  double total_density() const noexcept { return d[8]; }
  double mean_molecular_mass() const noexcept { return d[9]; }
  double &hydrogen() noexcept { return d[0]; }
  double &helium() noexcept { return d[1]; }
  double &atomic_oxygen() noexcept { return d[2]; }
  double &molecular_nitrogen() noexcept { return d[3]; }
  double &molecular_oxygen() noexcept { return d[4]; }
  double &atomic_nitrogen() noexcept { return d[5]; }
  double &temperature() noexcept { return d[6]; }
  double &exospheric_temperature() noexcept { return d[7]; }
  double &total_density() noexcept { return d[8]; }
  double &mean_molecular_mass() noexcept { return d[9]; }
}; // Dtm2020Output
} // namespace details

class Dtm2020 {
  static constexpr const int nlatm = 96;
  details::Dtm2020Input in;
  details::Dtm2020Output out;
  SolarActivityData flux_data;
  double hl0, ch, sh, c2h, s2h, c3h, s3h;
  double p10, p20, p30, p40, p50, p60, p11, p21, p31, p41, p51, p22, p32, p42,
      p52, p62, p33, p10mg, p20mg, p40mg, p11mg, p22mg, p31mg, p30mg, p50mg,
      p60mg;
  /* solar flux data (set by using get_solar_data) */
  double mf[2], mfbar[2], makp[4];
  /* data file values from DTM_2020_F107_Kp.dat */
  double az[nlatm], az2[nlatm], h[nlatm], he[nlatm], o[nlatm], o2[nlatm],
      t0[nlatm], tp[nlatm], tt[nlatm];

  int map_model_coeffs(const char *fn) noexcept;
  int get_solar_data(const dso::TwoPartDate &tutc,
                     bool use_adjusted_values = true) noexcept;

public:
  Dtm2020(const char *fn);

  void set_flux_data(const SolarActivityData &fdata) noexcept {flux_data=fdata;}

  /*
   * @param[in] longitude Longitude [-π, π] [rad]
   * @param[in] latitude  Latitude [-π/2, π/2] [rad]
   * @param[in] latitude  Altitude [km]. Should be > 120 km
   * @param[in] t         UTC datetime
   * @return Anything other than zero, means that the input parameters are
   *         erronuous
   */
  int set(const dso::TwoPartDate &t, double longitude, double latitude,
          double altitude_km) noexcept {
    in.set_position(longitude, latitude, altitude_km);
    in.set_date(t);
    return (in.validate_params() && get_solar_data(t));
  }

  double total_density_grcm3() const noexcept { return out.total_density(); }
  double total_density_kgm3() const noexcept {
    return out.total_density() * 1e-3 * 1e-4;
  }

  double temperature() const noexcept {return out.temperature();}
  double exospheric_temperature() const noexcept {return out.exospheric_temperature();}

  /* core function */
  int dtm3() noexcept;

  /* function heavily used to within the core of DTM2020 */
  double gldtm(const double *__restrict__ a, double *__restrict__ da,
               double ff0, double longitude) noexcept;

  int debug_set(double _lon, double _lat, double _alt, double _lh, double _doy, const double *f,
                const double *fbar, const double *akp) {
    in.lon = _lon;
    in.lat = _lat;
    in.alt = _alt;
    in.lh = _lh;
    in.doy = _doy;
    mf[0] = f[0]; mf[1] = f[1];
    mfbar[0] = fbar[0]; mfbar[1] = fbar[1];
    makp[0] = akp[0]; makp[1] = akp[1]; makp[2] = akp[2]; makp[3] = akp[3];
    return 0;
    }

}; // class Dtm2020

} // namespace dso
#endif
