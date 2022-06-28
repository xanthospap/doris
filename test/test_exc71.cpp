#include "astrodynamics.hpp"
#include "egravity.hpp"
#include "icgemio.hpp"
#include "geodesy/units.hpp"
#include "eigen3/Eigen/Eigen"
#include <matvec/vector3d.hpp>
#include <cassert>

constexpr const int Degree = 20;
constexpr const int Order = 20;

struct AuxParams {
  dso::HarmonicCoeffs *hc;
  dso::Mat2D<dso::MatrixStorageType::Trapezoid> *V, *W;
  int degree,order;
};

// handle gravity filed
int gravity(const char *gfn, int degree, int order, dso::HarmonicCoeffs &hc) {
  dso::Icgem gfc(gfn);
  // parse the header ...
  if (gfc.parse_header()) {
    fprintf(stderr, "ERROR! Failed to parse icgem header!\n");
    return 1;
  }

  assert(degree <= gfc.degree() && order <= degree && hc.degree() == degree);
  printf("Note: Setting degree for spherical harmonics to %d\n", degree);
  printf("Note: Setting order for spherical harmonics to %d\n",order);

  // parse data; store coefficients to hc
  if (gfc.parse_data(degree, order, &hc)) {
    fprintf(stderr, "ERROR! Failed to parse harmonic coefficients\n");
    return 1;
  }

  // de-normalize the harmonics coeffs
  hc.denormalize();

  return 0;
}

// Computes the variational equations
void variational_equations(double t,
                           // state and state transition matrix
                           const double *yNphi,
                           // state derivative and state transition matrix derivative
                           double *yNphi_derivs,
                           //
                           void *pAux
{

  printf("--- Variational Equations ---\n");
  assert(t==0);
  
  Eigen::Matrix<double, 3, 1> r; r << yNphi[0], yNphi[1], yNphi[2];
  printf("Satellite position  : %+15.6f %+15.6f %+15.6f\n", r(0), r(1), r(2));
  Eigen::Matrix<double, 3, 1> v; v << yNphi[3], yNphi[4], yNphi[5];
  printf("Satellite velocity  : %+15.6f %+15.6f %+15.6f\n", v(0), v(1), v(2));

  // void pointer to AuxParameters
  AuxParams *params = static_cast<AuxParams*>(pAux);

  // compute gravity-induced acceleration and gradient
  Eigen::Matrix<double, 3, 3> gpartials;
  Eigen::Matrix<double, 3, 1> gacc =
      dso::grav_potential_accel(r, params->degree, params->order, *(params->V),
                                *(params->W), *(params->hc), gpartials);

  printf("Gravity Acceleration: %+15.6f %+15.6f %+15.6f\n", gacc(0), gacc(1), gacc(2));
  printf("Gravity Gradient    : %+15.6f %+15.6f %+15.6f\n", gpartials(0,0), gpartials(0,1), gpartials(0,2));
  printf("                    : %+15.6f %+15.6f %+15.6f\n", gpartials(1,0), gpartials(1,1), gpartials(1,2));
  printf("                    : %+15.6f %+15.6f %+15.6f\n", gpartials(2,0), gpartials(2,1), gpartials(2,2));
  
  // state derivative (aka [v,a])
  Eigen::Matrix<double, 6, 1> dy;
  dy.block<3,1>(0,0) = v;
  dy.block<3,1>(3,0) = gacc;

  // derivative of state transition matrix
  Eigen::Matrix<double,6,6> dfdy;
  dfdy.block<3,3>(0,0) = Eigen::Matrix<double,3,3>::Zero();
  dfdy.block<3,3>(0,3) = Eigen::Matrix<double,3,3>::Identity();
  dfdy.block<3,3>(3,0) = gpartials;
  dfdy.block<3,3>(3,3) = Eigen::Matrix<double,3,3>::Zero();
  
  dPhi = dfdy * Phi;
  printf("Derivative of state transition matrix:\n");
  for (int i=0;i<6;i++) {
    for (int j=0;j<6;j++) {
      printf("%+15.6f  ", dPhi(i,j));
    }
    printf("\n");
  }

  // Derivative of combined state vector and state transition matrix
}


// Computes the variational equations
void variational_equations(double t,
                           // state
                           const Eigen::Matrix<double, 6, 1> &y,
                           // state transition matrix
                           const Eigen::Matrix<double, 6, 6> &Phi,
                           const AuxParams &params,
                           // state derivative
                           Eigen::Matrix<double, 6, 1> &dy,
                           // state transition matrix derivative
                           Eigen::Matrix<double, 6, 6> &dPhi) {
  printf("--- Variational Equations ---\n");
  assert(t==0);
  
  Eigen::Matrix<double, 3, 1> r = y.block<3,1>(0,0);
  printf("Satellite position  : %+15.6f %+15.6f %+15.6f\n", r(0), r(1), r(2));
  Eigen::Matrix<double, 3, 1> v = y.block<3,1>(3,0);
  printf("Satellite velocity  : %+15.6f %+15.6f %+15.6f\n", v(0), v(1), v(2));

  // compute gravity-induced acceleration and gradient
  Eigen::Matrix<double, 3, 3> gpartials;
  Eigen::Matrix<double, 3, 1> gacc =
      dso::grav_potential_accel(r, params.degree, params.order, *(params.V),
                                *(params.W), *(params.hc), gpartials);

  // state derivative (aka [v,a])
  dy.block<3,1>(0,0) = v;
  dy.block<3,1>(3,0) = gacc;
  
  assert(r(0)==y(0) && r(1)==y(1) && r(2) == y(2));
  printf("Gravity Acceleration: %+15.6f %+15.6f %+15.6f\n", gacc(0), gacc(1), gacc(2));
  printf("Gravity Gradient    : %+15.6f %+15.6f %+15.6f\n", gpartials(0,0), gpartials(0,1), gpartials(0,2));
  printf("                    : %+15.6f %+15.6f %+15.6f\n", gpartials(1,0), gpartials(1,1), gpartials(1,2));
  printf("                    : %+15.6f %+15.6f %+15.6f\n", gpartials(2,0), gpartials(2,1), gpartials(2,2));

  // derivative of state transition matrix
  Eigen::Matrix<double,6,6> dfdy;
  dfdy.block<3,3>(0,0) = Eigen::Matrix<double,3,3>::Zero();
  dfdy.block<3,3>(0,3) = Eigen::Matrix<double,3,3>::Identity();
  dfdy.block<3,3>(3,0) = gpartials;
  dfdy.block<3,3>(3,3) = Eigen::Matrix<double,3,3>::Zero();

  dPhi = dfdy * Phi;
  printf("Derivative of state transition matrix:\n");
  for (int i=0;i<6;i++) {
    for (int j=0;j<6;j++) {
      printf("%+15.6f  ", dPhi(i,j));
    }
    printf("\n");
  }

  // Derivative of combined state vector and state transition matrix

}

int main(int argc, char *argv[]) {
  if (argc < 2 || argc > 4) {
    fprintf(stderr,
            "Usage: %s <GRAVITY MODEL FILE> [DEGREE - optional] [ORDER - "
            "optional]\n",
            argv[0]);
    return 1;
  }
  /* parse degree */
  int degree = Degree;
  if (argc >= 3)
    degree = std::atoi(argv[2]);
  /* parse order */
  int order = Order;
  if (argc == 4)
    order = std::atoi(argv[3]);

  // Harmonic coefficients
  dso::HarmonicCoeffs hc(degree);
  if (gravity(argv[1], degree, order, hc)) {
    fprintf(stderr, "Aborting ...\n");
    return 1;
  }

  // Lagrange polynomials (depend on position) N+1 (for zero offset) and +2
  // because we are computing potential partials
  dso::Mat2D<dso::MatrixStorageType::Trapezoid> V(degree + 3, order + 3),
      W(degree + 3, order + 3);

  // Create Auxiliary parameters
  AuxParams params;
  params.hc = &hc;
  params.V = &V;
  params.W = &W;
  params.degree = degree;
  params.order = order;
  
  dso::OrbitalElements ele;
  ele.semimajor() = 6378.137e3 + 650.0e3;
  ele.eccentricity() = 0.001;
  ele.inclination() = dso::deg2rad(51e0);

  int error;
  auto y0 = dso::elements2state(hc.GM(), 0e0, ele, error);
  if (error) {
    fprintf(stderr, "ERROR. Failed to transform elements to state vector\n");
    return 1;
  }
  printf("State: ");
  for (int i=0;i<6;i++) printf("%+15.6f  ", y0(i));
  printf("\n");

  Eigen::Matrix<double, 7, 6> yPhi1;
  Eigen::Matrix<double, 6, 6> Phi = Eigen::Matrix<double, 6, 6>::Identity();
  
  yPhi1.block<1,6>(0,0) = y0; /* first row */
  yPhi1.block<6,6>(1,0) = Phi;
  for (int i=0; i<7; i++) {
    printf("  [%2d-%2d] %+15.3f %+15.3f %+15.3f %+15.3f %+15.3f %+15.3f\n",
    i*6, i*6+5, yPhi1(i,0), yPhi1(i,1), yPhi1(i,2), yPhi1(i,3), yPhi1(i,4), yPhi1(i,5));
  }

  // const auto yPhi2 = yPhi1;
  Eigen::Matrix<double, 6, 6> dPhi;
  Eigen::Matrix<double, 6, 1> dy;
  printf("State: ");
  for (int i=0;i<6;i++) printf("%+15.6f  ", y0(i));
  printf("\n");
  variational_equations(0e0,y0,Phi,params,dy,dPhi);

  return 0;
}
