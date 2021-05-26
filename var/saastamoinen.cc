/// @brief Gravity correction function for Saastamoinen model
/// @param[in] lat The geodetic latitude of the station (radians)
/// @param[in] H   the geodetic height of the station in meters
/// @return Gravity correction term aka fs(φ,H)
/// @see IERS conventions 2010, ch. 9.1.1, eq. 9.4 
inline double fs_phiH(double lat, double H) noexcept {
  return 1e0 − 0.00266e0*std::cos(2e0*lat) - 0.00000028e0 * H;
}

/// @brief Compute zenith hydrostatic delay in meters using Saastamoinen model
/// @param[in] lat The geodetic latitude of the station (radians)
/// @param[in] H   the geodetic height of the station in meters
/// @param[in] P0 total atmospheric pressure in hPa at the antenna reference 
///               point (e.g. antenna phase center for GPS)
/// @return The zenith hydrostatic delay in meters
/// @see IERS conventions 2010, ch. 9.2, eq. 9.11
double saastamoinen_hz(double lat, double H, double P0) noexcept {
  return (0.0022768e0 * P0) / fs_phiH(lat, H);
}

/// @brief Compute zenith wet delay in meters using Saastamoinen model
/// @param[in] e0 Partial pressure of water vapour at the ground beacon (mbar 
///               or hPa)
/// @param[in] T0 Surface temperature in Kelvin
/// @return wet zenith delay in meters
double saastamoinen_wet(double e0, double T0) noexcept {
  return 0.0022768e0 * (1255e0/T0 + 0.053e0)*e0;
}
