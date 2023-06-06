#include "astrodynamics.hpp"

/*
 * @param[in] sv The satellite's macromodel (components), where for each
 *               surface the normal vector is w.r.t the ECI frame, NOT the
 *               SV-fixed frame
 * @param[in] r  Satellite position vector in ECI frame
 * @param[in] sun Sun vector in ECI frame
 * @param[in] Mass Satellite's mass
 */
Eigen::Matrix<double, 3, 1>
dso::direct_solar_radiation_pressure(const std::vector<dso::MacroModelComponent> &sv,
                                const Eigen::Matrix<double, 3, 1> &r,
                                const Eigen::Matrix<double, 3, 1> &sun,
                                double Mass) noexcept {
  Eigen::Matrix<double, 3, 1> S = Eigen::Matrix<double, 3, 1>::Zero();
  
  /* Sun to satellite unit vector */
  const Eigen::Matrix<double, 3, 1> s = (sun - r).normalized();
  for (const auto &plate : sv) {
    /* surface normal, w.r.t ECI */
    const Eigen::Matrix<double, 3, 1> n = plate.normal().normalized();
    /* cosθ = <n,sun> */
    const double ct = n.dot(s);
    Eigen::Matrix<double, 3, 1> Ai = Eigen::Matrix<double, 3, 1>::Zero();
    /* specular */
    Ai += 2e0 * plate.optical_specular() * ct * n;
    /* diffusion */
    Ai += plate.optical_diffusion() * (s - (2e0 / 3e0) * n);
    /* absorption */
    Ai += plate.optical_absorption() * s;
    /* visibility */
    const double visibility = (ct >= 0e0) * 1e0;
    S += (Ai * plate.area() * ct) * visibility;
  }

  /* scale with (1/M) * (W/c) */
  constexpr const double P = 4.56316e-6;
  constexpr const double AU = iers2010::AU;
  const double R = (sun - r).norm();
  S *= (1e0 / Mass) * P * AU * AU / R;

  return S;
}
