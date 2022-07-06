#include <cstdio>
#include "sp3/sp3.hpp"
#include "astrodynamics.hpp"
#include "egravity.hpp"
#include "icgemio.hpp"
#include "geodesy/units.hpp"
#include "eigen3/Eigen/Eigen"
#include "integrators.hpp"
#include "datetime/datetime_write.hpp"

using dso::sp3::SatelliteId;
constexpr const int Degree = 20;
constexpr const int Order = 20;

// to transfer parameters for Variational Equations
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

  // parse data; store coefficients to hc
  if (gfc.parse_data(degree, order, &hc)) {
    fprintf(stderr, "ERROR! Failed to parse harmonic coefficients\n");
    return 1;
  }

  // de-normalize the harmonics coeffs
  hc.denormalize();

  return 0;
}

void VariationalEquations(
    [[maybe_unused]] double t,
    // state and state transition matrix
    const Eigen::VectorXd &yPhi,
    // state derivative and state transition matrix derivative
    Eigen::Ref<Eigen::VectorXd> yPhiP,
    //
    void *pAux) noexcept 
{
  static unsigned call_nr = 0;

  // void pointer to AuxParameters
  AuxParams *params = static_cast<AuxParams *>(pAux);

  // split position and velocity vectors
  Eigen::Matrix<double, 3, 1> r = yPhi.block<3, 1>(0, 0);
  Eigen::Matrix<double, 3, 1> v = yPhi.block<3, 1>(3, 0);

  // compute gravity-induced acceleration and gradient
  Eigen::Matrix<double, 3, 3> gpartials;
  Eigen::Matrix<double, 3, 1> gacc =
      dso::grav_potential_accel(r, params->degree, params->order, *(params->V),
                                *(params->W), *(params->hc), gpartials);

  // State transition (skip first column which is the state vector)
  Eigen::Matrix<double,6,6> Phi (yPhi.data()+6);

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

  // matrix to vector (column-wise)
  yPhiP = Eigen::VectorXd(
      Eigen::Map<Eigen::VectorXd>(yPhip.data(), yPhip.cols() * yPhip.rows()));

  ++call_nr;
  return;
}

int main(int argc, char *argv[]) {

  // handle gravity field
  if (argc < 3 || argc > 5) {
    fprintf(stderr,
            "Usage: %s <SP3 FILE> <GRAVITY MODEL FILE> [DEGREE - optional] "
            "[ORDER - optional]\n",
            argv[0]);
    return 1;
  }

  // parse degree and order of gravity field
  int degree = Degree;
  int order  = Order;
  if (argc >= 4)
    degree = std::atoi(argv[3]);
  if (argc == 5)
    order = std::atoi(argv[4]);

  // Harmonic coefficients
  dso::HarmonicCoeffs hc(degree);
  if (gravity(argv[2], degree, order, hc)) {
    fprintf(stderr,
            "ERROR. Failed to compute harmonic coefficients. Aborting ...\n");
    return 1;
  }

  // Lagrange polynomials (depend on position) N+1 (for zero offset) and +2
  // because we are computing potential partials. Allocate space ...
  dso::Mat2D<dso::MatrixStorageType::Trapezoid> V(degree + 3, order + 3),
      W(degree + 3, order + 3);

  // Handle the Sp3 file, prepare for parsing ...
  dso::Sp3c sp3(argv[1]);
  #ifdef DEBUG
  sp3.print_members();
  #endif

  SatelliteId sv;
  dso::Sp3DataBlock block;

  if (sp3.num_sats() == 1) {
    sv.set_id(sp3.sattellite_vector()[0].id);
  } else {
    fprintf(stderr, "More than one Satellites in sp3 file. This test program "
                    "is only meant to work with one.\n");
    return 1;
  }

  // Create Auxiliary parameters
  AuxParams params{&hc, &V, &W, degree, order};
  // SetUp an integrator
  const double relerr = 1.0e-13; // Relative and absolute
  const double abserr = 1.0e-6;  // accuracy requirement
  dso::SGOde Integrator(VariationalEquations, /*neqn*/ 6 * 7, relerr, abserr,
                        &params);
  Eigen::Matrix<double,6*7,1> yPhi; 
  Eigen::VectorXd sol(6*7);
  Eigen::VectorXd r0(6);
  // Time
  const auto t0 = sp3.start_epoch();
  auto epoch = t0;
  const auto step = sp3.interval();

  // let's try reading the records; note that -1 denotes EOF
  int error = 0;
  std::size_t rec_count = 0;
  char buf[64];
  do {
    
    // read nex record ...
    if ((error=sp3.get_next_data_block(sv, block))>0) {
      printf("Something went wrong ....status = %3d\n", error);
      return 1;
    } 
  
    // check the heath status
    bool position_ok = !block.flag.is_set(dso::Sp3Event::bad_abscent_position);
    
    if (position_ok) {
      yPhi = Eigen::Matrix<double,6*7,1>::Zero();

      // accumulate state (m, m/sec)
      r0 << block.state[0]*1e3, block.state[1]*1e3, block.state[2]*1e3,
        block.state[4]*1e-2, block.state[5]*1e-2, block.state[6]*1e-2;
      yPhi.block<6,1>(0,0) = r0;
      
      // set the Phi part (aka, the state transition matrix) to I(6x6)
      int offset = 0;
      for (int col=1; col<7; col++) {
        int index = 6*col + offset++;
        yPhi(index) = 1e0;
      }

      // t = current_time - t0 [seconds]
      double t = block.t.delta_sec(t0).to_fractional_seconds();
      // tout = t + step [seconds]
      double tout = t + step.to_fractional_seconds();
      
      // integrate
      Integrator.flag() = 1;
      Integrator.de(t, tout, yPhi, sol);

      // output epoch as datetime
      epoch = block.t;
      epoch.add_seconds( dso::milliseconds(static_cast<long>(tout*1e3)) );
      dso::strftime_ymd_hmfs(epoch, buf);

      printf("\n%s %+15.4f %+15.4f %+15.4f %+15.7f %+15.7f %+15.7f\n", 
        buf, sol(0), sol(1), sol(2), sol(3), sol(4), sol(5));
        
    }
    ++rec_count;
  
  } while (!error);

  printf("Num of records read: %6lu\n", rec_count);

  return (error==-1);
}
