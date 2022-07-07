#include "astrodynamics.hpp"
#include "egravity.hpp"
#include "icgemio.hpp"
#include "geodesy/units.hpp"
#include "eigen3/Eigen/Eigen"
#include <matvec/vector3d.hpp>
#include <cassert>
#include "integrators.hpp"

constexpr const int Degree = 20;
constexpr const int Order = 20;
constexpr const double GM = 398600.4415e+9;
// constexpr const double Re = 6378.137e3;

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
void variational_equations(
    [[maybe_unused]] double t,
    // state and state transition matrix
    const Eigen::VectorXd &yPhi,
    // state derivative and state transition matrix derivative
    Eigen::Ref<Eigen::VectorXd> yPhiP,
    //
    void *pAux) noexcept {

  static unsigned call_nr = 0;
  // printf("--- Variational Equations --- (%u)\n", call_nr);
  //assert(t == 0);

  Eigen::Matrix<double, 3, 1> r = yPhi.block<3, 1>(0, 0);
  Eigen::Matrix<double, 3, 1> v = yPhi.block<3, 1>(3, 0);
  //printf("%u ", call_nr);
  //for (int i=0; i<3; i++) printf("%+20.10f ", r(i));
  //for (int i=0; i<3; i++) printf("%+20.10f ", v(i));
  //printf("\n");

  // void pointer to AuxParameters
  AuxParams *params = static_cast<AuxParams *>(pAux);

  // compute gravity-induced acceleration and gradient
  Eigen::Matrix<double, 3, 3> gpartials;
  // for consistency with Montenbrck, use
  Eigen::Matrix<double, 3, 1> gacc = dso::grav_potential_accel(
      r, params->degree, params->order, *(params->V),
      *(params->W), *(params->hc), gpartials);

  //printf("%+15.10f %+15.10f %+15.10f\n", gacc(0), gacc(1),
  //       gacc(2));
  //printf("Gravity Gradient    : %+15.6f %+15.6f %+15.6f\n", gpartials(0, 0),
  //       gpartials(0, 1), gpartials(0, 2));
  //printf("                    : %+15.6f %+15.6f %+15.6f\n", gpartials(1, 0),
  //       gpartials(1, 1), gpartials(1, 2));
  //printf("                    : %+15.6f %+15.6f %+15.6f\n", gpartials(2, 0),
  //       gpartials(2, 1), gpartials(2, 2));

  // State transition (skip first column which is the state vector)
  Eigen::Matrix<double,6,6> Phi (yPhi.data()+6);
  //for (int i=0; i<36; i++) {
  //  printf("%20.15f ", Phi.data()[i]);
  //}
  //printf("\n");

  // derivative of state transition matrix, aka
  // | v (3x3)   0 (3x3)     I (3x3)   |
  // | a (3x1) da/dr (3x3) da/dv (3x3) |
  // because dv/dr = 0 and
  //         da/dr = I
  Eigen::Matrix<double,6,6> dfdy;
  dfdy.block<3, 3>(0, 0) = Eigen::Matrix<double, 3, 3>::Zero();
  dfdy.block<3, 3>(0, 3) = Eigen::Matrix<double, 3, 3>::Identity();
  dfdy.block<3, 3>(3, 0) = gpartials;
  dfdy.block<3, 3>(3, 3) = Eigen::Matrix<double, 3, 3>::Zero();

  // Derivative of combined state vector and state transition matrix
  // dPhi = dfdy * Phi;
  Eigen::Matrix<double,6,7> yPhip;
  yPhip.block<6,6>(0,1) = dfdy * Phi;

  // state derivative (aka [v,a]), in one (first) column
  yPhip.block<3,1>(0,0) = v;
  yPhip.block<3,1>(3,0) = gacc;

  // yPhiP = Eigen::VectorXd(yPhip.data());
  yPhiP = Eigen::VectorXd(
      Eigen::Map<Eigen::VectorXd>(yPhip.data(), yPhip.cols() * yPhip.rows()));
 
  ++call_nr; 
  return;
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
  auto y0 = dso::elements2state(GM, 0e0, ele, error);
  if (error) {
    fprintf(stderr, "ERROR. Failed to transform elements to state vector\n");
    return 1;
  }
  //printf("State at t0 is: ");
  //for (int i=0;i<6;i++) printf("%+20.10f  ", y0(i));
  //printf("\n");

  Eigen::Matrix<double, 6,7> yPhiM;
  // first column is the state vector
  yPhiM.block<6,1>(0,0) = y0;
  yPhiM.block<6,6>(0,1) = Eigen::Matrix<double, 6, 6>::Identity();
  
  // arrange matrix to one big vector
  Eigen::Matrix<double, 6*7, 1> yPhi (yPhiM.data());
  
  // SetUp an integrator
  const double relerr = 1.0e-13;                     // Relative and absolute
  const double abserr = 1.0e-6;                      // accuracy requirement
  const double t0 = 0e0;
  const double t_end = 300.0; //86400.0;
  const double step = 300.0;

  dso::SGOde Sg(variational_equations, /*neqn*/6*7, relerr, abserr, &params);

  double t = t0;
  Eigen::VectorXd sol(6*7);
  while (t<t_end) {

    // output time
    double tout = t + step;

    // integrate 
    Sg.de(t, tout, yPhi, sol);

    // print solution
    printf("Solution at t=%.3f\n", t);
    printf("%+15.6f %+15.6f %+15.6f %+15.6f %+15.6f %+15.6f\n",
      sol(0),sol(1),sol(2),sol(3),sol(4),sol(5));

  }

  return 0;
}
