#ifndef __DSO_SAT_ATTITUDE_MODELS_HPP__
#define __DSO_SAT_ATTITUDE_MODELS_HPP__

#include "eigen3/Eigen/Eigen"
#include "satellites.hpp"
#include "satellites/jason3_quaternions.hpp"
#include "iers2010/rotations.hpp"
#include <cassert>

namespace dso {

class SvFrame {
  /* CoG point in SV-fixed RF */
  Eigen::Matrix<double, 3, 1> cog_sf;
  /* Antenna RP point in SV-fixed RF */
  Eigen::Matrix<double, 3, 1> arp_sf;
  /* quaternion hunter for attitude (body-frame) */
  dso::JasonQuaternionHunter qhunt;
  /* left/right array angles for attitude */
  dso::JasonSolarArrayHunter shunt;
  /* satellite macromodel */
  dso::MacroModel<SATELLITE::Jason3> mplates;
  double sat_mass{0e0};

public:
  SvFrame(Eigen::Matrix<double, 3, 1> vcog, Eigen::Matrix<double, 3, 1> varp,
          const char *qfn, const char *sfn, double mass)
      : cog_sf(vcog), arp_sf(varp), qhunt(qfn), shunt(sfn), sat_mass(mass) {}

  double mass() const noexcept {return sat_mass;}

  const dso::MacroModel<SATELLITE::Jason3> &plates() const noexcept {
    return mplates;
  }

  dso::iStatus macromodel_iterator(const dso::TwoPartDate &tai, std::vector<dso::MacroModelComponent> &sv) noexcept {
    sv.clear();
    sv.reserve(mplates.NumPlates);

    /* get rotation quaternion for body frame */
    Eigen::Quaternion<double> q;
    if (qhunt.get_at(tai, q)) {
      fprintf(
          stderr,
          "[ERROR] Failed to find quaternion for datetime (traceback: %s)\n",
          __func__);
      return dso::iStatus(1);
    }
    /* get rotation for left and right arrays */
    double lrot, rrot;
    if (shunt.get_at(tai, lrot, rrot)) {
      fprintf(stderr,
              "[ERROR] Failed getting solar array rotation (traceback: %s)\n",
              __func__);
      return dso::iStatus(1);
    }
    /* iterate macromodel; for sv body */
    for (int i=0; i<mplates.NumPlates-2; i++) {
      dso::MacroModelComponent plate(mplates.plate(i));
      /* body-fixed to ECI */
      plate.set_normal(q.conjugate().normalized() * mplates.plate(i).normal());
      /* push back */
      sv.emplace_back(plate);
    }

    /* left solar array */
    {
      const int index = 6;
      dso::MacroModelComponent plate(mplates.plate(index));
      /* left solar array default normal */
      Eigen::Matrix<double,3,1> n = mplates.plate(index).normal();
      /* rotate solar array to body-fixed frame */
      n = dso::rotate<dso::RotationAxis::Y>(lrot, n);
      /* body-fixed to ECI */
      plate.set_normal(q.conjugate().normalized() * n);
      /* push back */
      sv.emplace_back(plate);
    }
    
    /* right solar array */
    {
      const int index = 7;
      dso::MacroModelComponent plate(mplates.plate(index));
      /* right solar array default normal */
      Eigen::Matrix<double,3,1> n = mplates.plate(index).normal();
      /* rotate solar array to body-fixed frame */
      n = dso::rotate<dso::RotationAxis::Y>(rrot, n);
      /* body-fixed to ECI */
      plate.set_normal(q.conjugate().normalized() * n);
      /* push back */
      sv.emplace_back(plate);
    }

    return dso::iStatus(0);
  }

  /* body-fixed frame to GCRF */
  int arp2cog(const dso::TwoPartDate &tai, Eigen::Matrix<double, 3, 1> &dr) {
    /* r_sat_fixed = Q * r_sat_gcrs */
    Eigen::Quaternion<double> q;
    if (qhunt.get_at(tai, q)) {
      fprintf(stderr, "[ERROR] Failed to find quaternion for datetime\n");
      return 1;
    }
    dr = q.conjugate().normalized() * (cog_sf - arp_sf);
    return 0;
  }

  /*
  int get_attitude_quaternion(const dso::TwoPartDate &tai, Eigen::Quaternion<double> &q) noexcept {
    if (qhunt.get_at(tai, q)) {
      fprintf(stderr, "[ERROR] Failed to find quaternion for datetime\n");
      return 1;
    }
    return 0;
  }
  */

  ///* returns (S / m) where S is the projected area */
  //int projected_area(const dso::TwoPartDate &tai,
  //                   const Eigen::Matrix<double, 3, 1> &vrel, double &b) noexcept {
  //  /* attitude matrix */
  //  Eigen::Quaternion<double> q;
  //  if (qhunt.get_at(tai, q)) {
  //    fprintf(stderr, "[ERROR] Failed to find quaternion for datetime\n");
  //    return 1;
  //  }
  //  /* normalized relative velocity */
  //  const auto v = vrel.normalized();
  //  /* iterate through SV plates */
  //  double S = 0e0;
  //  const int numPlates = mplates.NumPlates;
  //  const MacroModelComponent *plate = nullptr;
  //  for (int i=0; i<numPlates; i++) {
  //    plate = mplates.mmcomponents + i;
  //    /* (unit) vector of plate, normal in sv-fixed RF */
  //    Eigen::Matrix<double,3,1> n(plate->m_normal);
  //    const double cosA = (q.conjugate().normalized()*n).transpose() * v;
  //    if (cosA > 0e0) {
  //      S += (plate->m_surf * cosA);
  //    }
  //  }
  //  /* compute ballistic coefficient, without the Cd */
  //  b = S/sat_mass;
  //  return 0;
  //}
}; /* SvFrame */

}/* namespace dso */

#endif
