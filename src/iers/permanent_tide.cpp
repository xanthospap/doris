namespace {
    /// @brief Computes the Lagrange function of order = 2, aka P_2, as
    ///        P2(sinφ) = (3 * sin^2(φ)-1) / 2
    /// @return 
    inline double p2f(double sinlat) noexcept {
        return (3e0 *sinlat * sinlat - 1e0) / 2e0;
    }    
}// unnamed namespace

/// @brief Compute site displacements to get “mean tide” coordinates from
///         “conventional tide free” coordinates. See Section 7.1.1.2 Permanent
///         deformation, of IERS 2010
/// @param sites List of sites, in geodetic ECEF coordinates [λ,φ,h] in [rad] 
///              and[m]
/// @param Senu  Corresponding corrections to be added to the “conventional tide
///              free” computed position to obtain the “mean tide” position. 
///              This vector contains the radial and traverse (northwards) 
///              elements as: [0,D_traverse,D_radial] in [m]
/// @return Always 0
int permanent_tide(const std::vector<Eigen::Matrix<double, 3, 1>> &sites,
                   std::vector<Eigen::Matrix<double, 3, 1>> &Senu) noexcept 
{
  if (Senu.capacity() < sites.size())
    Senu.reserve(sites.size());

  Senu.clear();

    for (const auto &site: sites) {
        const double lat = site(1);
        const double slat = std::sin(lat);
        const double s2lat = std::sin(2e0*lat);
        const double p2 = p2f(slat);
        // apply equations 14b and 14a
        Senu.push_back(
            Eigen::Matrix<double, 3, 1>(0e0, 
                (-0.0252e0 - 0.0001e0*p2) * s2lat), 
                (-0.1206e0 + 0.0001e0*p2)*p2);
    }

    return 0;
}
