#ifndef __DSO_NRLMSISE00_CVERS_HPP__
#define __DSO_NRLMSISE00_CVERS_HPP__

#include "datetime/dtcalendar.hpp"
#include "nrlmsise00_const.hpp"
#include "var_utils.hpp"
#include <cmath>
#include <cstdint>
#include <cstring>
#ifdef DEBUG
#include <cassert>
#endif

namespace dso {

namespace nrlmsise00 {

namespace detail {

/// @brief Used within NRLMSISE00
double ccor2(double alt, double r, double h1, double zh, double h2) noexcept;

/// @brief Used within NRLMSISE00
double ccor(double alt, double r, double h1, double zh) noexcept;

/// @brief Used within NRLMSISE00
double dnet(double &dd, double dm, double zhm, double xmm, double xm) noexcept;

/// @brief Used within NRLMSISE00
inline double g0(double a, const double *p) noexcept {
  return (a - 4e0 +
          (p[25] - 1e0) * (a - 4e0 +
                           (std::exp(-std::abs(p[24]) * (a - 4e0)) - 1e0) /
                               std::abs(p[24])));
}

/// @brief Used within NRLMSISE00
inline double sumex(double ex) noexcept {
  return 1e0 + (1e0 - std::pow(ex, 19e0)) / (1e0 - ex) * std::pow(ex, 0.5e0);
}

/// @brief Used within NRLMSISE00
inline double sg0(double ex, const double *p, const double *ap) noexcept {
  return (g0(ap[1], p) + (g0(ap[2], p) * ex + g0(ap[3], p) * ex * ex +
                          g0(ap[4], p) * std::pow(ex, 3e0) +
                          (g0(ap[5], p) * std::pow(ex, 4e0) +
                           g0(ap[6], p) * std::pow(ex, 12e0)) *
                              (1e0 - std::pow(ex, 8e0)) / (1e0 - ex))) /
         sumex(ex);
}

/// @brief Spline interpolation, used within NRLMSISE00
void spline(const double *__restrict__ x, const double *__restrict__ y, int n,
            double yp1, double ypn, double *__restrict__ y2,
            double *work /*size >= n */) noexcept;

/// @brief Spline interpolation, used within NRLMSISE00
double splini(const double *__restrict__ xa, const double *__restrict__ ya,
              double *__restrict__ y2a, int n, double x) noexcept;

/// @brief Spline interpolation, used within NRLMSISE00
double splint(const double *__restrict__ xa, const double *__restrict__ ya,
              double *__restrict__ y2a, int n, double x) noexcept;

/// @brief A struct to hold the input switches for the NRLMSISE00 model
struct Switches {
  static constexpr const int dim = 25;
  double sw[dim] = {0e0};
  double swc[dim] = {0e0};

  /// @brief Constructor sets all switches to on (aka=1)
  Switches() noexcept { set_on(); }

  /// @brief Set value of given switch; this will perform the following:
  /// sw[index] = val % 2 , and
  /// swc[index] = (val == 1 || val == 2) ? 1e0 : 0e0;
  /// As performed in the original FORTRAN code using TSELEC
  /// @param[in] index The switch index
  /// @param[in] val   Val to set (kinda)
  void set_switch(int index, int val) noexcept {
#ifdef DEBUG
    assert(index >= 0 && index < dim);
#endif
    sw[index] = val % 2;
    swc[index] = (val == 1 || val == 2) ? 1e0 : 0e0;
  }

  /// @brief Set switch/switches to zero.
  /// @param[in] index If < 0, sets all switches (sw and swc) to 0. If >= 0,
  ///            only sets sw[index] and swc[index] to zero
  void set_null(int index = -1) noexcept {
    if (index < 0) {
      // std::memset(isw, 0, sizeof(int8_t) * dim);
      std::memset(sw, 0, sizeof(double) * dim);
      std::memset(swc, 0, sizeof(double) * dim);
    } else {
#ifdef DEBUG
      assert(index < dim);
#endif
      // isw[index] = 0;
      sw[index] = 0e0;
      swc[index] = 0e0;
    }
  }

  /// @brief Set switch/switches to one.
  /// @param[in] index If < 0, sets all switches (sw and swc) to 1. If >= 0,
  ///            only sets sw[index] and swc[index] to one.
  void set_on(int index = -1) noexcept {
    if (index < 0) {
      // for (int i=0; i<dim; i++) isw[i] = 1;
      for (int i = 0; i < dim; i++)
        sw[i] = 1e0;
      for (int i = 0; i < dim; i++)
        swc[i] = 1e0;
    } else {
#ifdef DEBUG
      assert(index < dim);
#endif
      // isw[index] = 1;
      sw[index] = 1e0;
      swc[index] = 1e0;
    }
  }
}; // dso::nrlmsise00::Switches

/// @brief Holds Ap indexes for NRLMSISE model. These are:
/// (0) DAILY AP
/// (1) 3 HR AP INDEX FOR CURRENT TIME
/// (2) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
/// (3) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
/// (4) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
/// (5) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
///     TO CURRENT TIME
/// (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR TO CURRENT
///     TIME struct
/// This struct is supposed to be a member of the more general
/// dso::NrlMsise00Input class.
struct ApArray {
  double a[7];
  ApArray& operator=(const ApArray& other) noexcept {
    if (this != &other) {
      std::memcpy(a, other.a, sizeof(double) * 7);
    }
    return *this;
  }
  ApArray& operator=(const double *other) noexcept {
    if (this->a != other) {
      std::memcpy(a, other, sizeof(double) * 7);
    }
    return *this;
  }
}; // dso::nrlmsise00::ApArray

/// @brief A class to hold input parameters for NRLMSISE00.
/// This is the core class, providing access to the data. In readl case, POD
/// scenarios, the data (flux, etc) need to change fase, hence we are later on
/// going to wrap this class with a data feed.
/// @note To use the aparr,  sw.sw[8]=-1 must hold. To enable this, use the
///       function: InParamsCore::use_aparray();
/// @note To request results in m^3 and kg/m^3, set meters to on
///       (InParamsCore::meters_on())
/// @note UT, Local Time, and Longitude are used independently in the
///       model and are not of equal importance for every situation.
///       For the most physically realistic calculation these three
///       variables should be consistent (STL=SEC/3600+GLONG/15).
///       The Equation of Time departures from the above formula
///       for apparent local time can be included if available but
///       are of minor importance.
/// @note F107 and F107A values used to generate the model correspond
///       to the 10.7 cm radio flux at the actual distance of the Earth from
///       the Sun rather than the radio flux at 1 AU.
struct InParamsCore {
  nrlmsise00::detail::Switches sw;
  int year;     // year, currently ignored
  int doy;      // day of year
  double sec;   // seconds in day (UT)
  double alt;   // altitude in kilometers
  double glat;  // geodetic latitude [degrees]
  double glon;  // geodetic longitude  [degrees]
  double lst;   // local apparent solar time (hours), see note
  double f107A; // 81 day average of F10.7 flux (centered on doy)
  double f107;  // daily F10.7 flux for previous day
  double ap;    // magnetic index(daily)
  nrlmsise00::detail::ApArray aparr; // ap index array
  bool meters_{false};               // request result in m^3 and kg/m^3

  InParamsCore() noexcept {};
  InParamsCore(dso::modified_julian_day mjd, double secday) noexcept;

  /// @brief Trigger ise of AP index array
  void use_aparray() noexcept { sw.set_switch(8, -1); };

  /// @brief Request result in m^3 and kg/m^3
  bool meters() const noexcept { return meters_; }

  /// @brief Check if result in m^3 and kg/m^3 is requested
  void meters_on() noexcept { meters_ = true; }

  /// @brief Compute and return the fractional doy
  double fdoy() const noexcept { return (double)doy + sec / 86400e0; }

  /// @brief Set all witches to 0 (off) -- both sw and swc --
  void set_switches_off() noexcept { sw.set_null(); }

  /// @brief Set all switches to 1 (on) -- both sw and swc --
  void set_switches_on() noexcept { sw.set_on(); }

  /// @brief Set a given switch
  /// @param[in] index Switch index
  /// @param[in] val   Value to set; The actual value set, is dominated by
  ///            Switches::set_switch
  void set_switch(int index, double val) noexcept { sw.set_switch(index, val); }

  /// @brief Set the Ap index array values using an input array
  /// @param[in] ap_array Should hold 7 double values, as ordered in ApArray
  void set_ap_array(const double *ap_array) noexcept {
    aparr = ap_array;
  }

  /// @brief Set the Ap index array values using an ApArray
  void set_ap_array(const ApArray &ap_array) noexcept {
    aparr = ap_array;
  }

  /// @brief
  int get_flux_data(const char *csvfn) noexcept;
}; // InParams

/// @brief The Nrlmsise00DataFeed class is meant to be a data provider for the
///         InParamsCore class.
/// It is meant to act as an "intermidiate" between the InParamsCore class and
/// a ClesTrak, Space Weather CSV file. It will parse the file, and extract
/// flux and Ap information depending on given dates. It is helpful usually in
/// the case of POD, where we need to have very dense and forward-in-time data
/// values for computing drag force.
struct Nrlmsise00DataFeed {
  ///< ClesTrak, Space Weather CSV file
  const char *fncsv_;
  ///< Current hour index [0-8), 3-hour intervals
  int chi;
  dso::modified_julian_day mjd_;
  ///< Flux and ap data for current and 3 days prior to current dates. Cuurent
  ///< date is always flux_data_[3] (one prior is flux_data_[2], two days prior
  /// is flux_data_[1], etc).
  dso::utils::celestrak::details::CelestTrakSWFlux flux_data_[4];
  ///< last position in file for which data were retrieved
  std::ifstream::pos_type fpos;

  /// @brief Constructor from ClesTrak, Space Weather CSV file
  Nrlmsise00DataFeed(const char *fncsv, InParamsCore &p);

  /// @brief Initialize -- do not use, is used automatically --
  int init(InParamsCore &p) noexcept;

  /// @brief Update flux/Ap records of InParamsCore, using its date/time data
  /// To update the flux/Ap data in an InParamsCore instance, set the target
  /// data (aka year and doy) and time (via sec), and pass it in to this
  /// function.
  int update(InParamsCore &p) noexcept;
};

/// @brief Enumeration to hold the type of data feed for NRLMSISE00 input
///        parameters.
enum class FluxDataFeedType : char {
  NONE,     ///< Manual data feed, user is responsible
  ST_CSV_SW ///< Automatic data feed using a ClesTrak, Space Weather CSV file
};

} // namespace detail

/// @brief A class to hold NRLMSISE00 input parameters and the way they are
///        retrieved/updated
template <detail::FluxDataFeedType T> struct InParams {};

/// @brief Specialization; no data feed, user is responsible for retrieving and
///        updating flux/Ap records
template <> struct InParams<detail::FluxDataFeedType::NONE> {
  detail::InParamsCore params_;
};

/// @brief Specialization; flux/Ap data are handled by an
/// detail::Nrlmsise00DataFeed instance, which is responsible for retrieving
/// and updating data, base on a CSV file
template <> struct InParams<detail::FluxDataFeedType::ST_CSV_SW> {
  detail::InParamsCore params_;
  detail::Nrlmsise00DataFeed feed_;

  InParams(const char *csv, dso::modified_julian_day mjd, double secday)
      : params_(mjd, secday), feed_(csv, params_) {}

  /// @brief Update datetime parameters for call to NRLMSISE. year, doy, sec,
  ///        are set and flux/Ap data are extracted from the CSV file.
  int update_params(int mjd, double secday) noexcept {
    const dso::datetime<dso::milliseconds> t(
        mjd, dso::milliseconds(static_cast<dso::milliseconds::underlying_type>(
                 secday * 1e3)));
    const auto ydoy = t.as_ydoy();
    params_.year = ydoy.__year.as_underlying_type();
    params_.doy = ydoy.__doy.as_underlying_type();
    params_.sec = t.sec().to_fractional_seconds();
    return feed_.update(params_);
  }
};

/// @brief NRLMSISE00 output parameters
struct OutParams {
  double d[9]; ///< densities
  double t[2]; ///< temperatures
};             // OutParams

} // namespace nrlmsise00

/// @brief A class to handle calls to the NRLMSISE00 model
class Nrlmsise00 {
  static const double ptm[10];
  static const double pdm[8][10];
  static const double pavgm[10];

  // two dim arrays
  double pma[10][100];
  double pd[9][150];
  double ptl[4][100];
  double plg[4][9];
  const double pdl[2][25];

  // one dim arrays ...
  double pt[150];
  double ps[150];
  double tn1[5];
  double tn2[4];
  double tn3[5];
  double tgn1[2];
  double tgn2[2];
  double tgn3[2];
  double apt[4];
  const double zn2[4] = {72.5e0, 55e0, 45e0, 32.5e0};
  const double zn3[5] = {32.5e0, 20e0, 15e0, 10e0, 0e0};
  // const double sam[100]; -> just a function

  double gsurf;
  double re;
  double dd;
  double dm28, tlb, xg0, dm28m;
  double dfa;
  double ctloc, stloc;
  double c2tloc, s2tloc;
  double s3tloc, c3tloc;
  double apdf;

  /// @brief Replaces the original FORTRAN implementation SAM array
  constexpr double sam(int i) const noexcept {
#ifdef DEBUG
    assert(i >= 0 && i < 100);
#endif
    return (i < 50);
  }

  double densm(double alt, double d0, double xm, double &tz) const noexcept;
  double densu(double alt, double dlb, double t1, double t2, double xm,
               double xalpha, double &tz, double zlb, double s2,
               const double *zn1) noexcept;
  double zeta(double zz, double zl) const noexcept {
    return (zz - zl) * (re + zl) / (re + zz);
  }
  double scalh(double alt, double xm, double temp) const noexcept {
    const double g = gsurf / std::pow(1e0 + alt / re, 2e0);
    return dso::nrlmsise00::detail::r100gas * temp / (g * xm);
  }
  double glatf(double lat, double &gv) const noexcept;

  double glob7s(const nrlmsise00::detail::InParamsCore *in, double *p) noexcept;

  int ghp7(const nrlmsise00::detail::InParamsCore *in,
           nrlmsise00::OutParams *out, double press) noexcept;

  int gts7(const nrlmsise00::detail::InParamsCore *in,
           nrlmsise00::OutParams *out, int mass = 48) noexcept;

  double globe7(const nrlmsise00::detail::InParamsCore *in, double *p) noexcept;

public:
  Nrlmsise00() noexcept;

  /// @brief Neutral Atmosphere Empirical Model from the surface to lower
  /// exosphere.
  ///
  /// @param[in] in Includes:
  ///   * doy day of year
  ///   * sec seconds in day
  ///   * alt altitude above surface (km)
  ///   * g_lat geodetic latitude
  ///   * g_long geodetic longitude
  ///   * lst local apparent solar time (hours)
  ///   * f107A 81 day average of F10.7 flux (centered on doy)
  ///   * f107 daily F10.7 flux for previous day
  ///   * ap magnetic index array
  ///   * d density array
  ///   * t temperature array
  ///
  /// The output variables are:
  ///  - out.d[0]: HE number density(cm-3)
  ///  - out.d[1]: O  number density(cm-3)
  ///  - out.d[2]: N2 number density(cm-3)
  ///  - out.d[3]: O2 number density(cm-3)
  ///  - out.d[4]: AR number density(cm-3)
  ///  - out.d[5]: total mass density(gm/cm3) [includes d[8] in td7d]
  ///  - out.d[6]: H Number density(cm-3)
  ///  - out.d[7]: N Number density(cm-3)
  ///  - out.d[8]: Anomalous oxygen number density(cm-3)
  ///  - out.t[0]: exospheric temperature
  ///  - out.t[1]: temperature at alt
  ///
  /// @note d[5] is the sum of the mass densities of the
  ///        species labeled by indices 0-4 and 6-7 in output variable d.
  ///        This includes He, O, N2, O2, Ar, H, and N but does NOT include
  ///        anomalous oxygen (species index 8).
  int gtd7(const nrlmsise00::detail::InParamsCore *in,
           nrlmsise00::OutParams *out, int mass = 48) noexcept;

  /// @brief Neutral Atmosphere Empirical Model from the surface to lower
  /// exosphere, including the anomalous oxygen contribution.
  ///
  ///  This subroutine provides Effective Total Mass Density for output
  ///  d[5] which includes contributions from "anomalous oxygen" which can
  ///  affect satellite drag above 500 km.
  int gtd7d(const nrlmsise00::detail::InParamsCore *in,
            nrlmsise00::OutParams *out, int mass = 48) noexcept;

}; // Nrlmsise00

} // namespace dso
#endif
