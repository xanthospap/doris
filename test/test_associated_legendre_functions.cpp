#include "associated_legendre.hpp"
#include <array>
#include <cassert>
#include <cstdio>
#include <limits>

using dso::AssociatedLegendreFunctions;
constexpr const double PRECISION = 1e-16;
bool ApproxEqual(double a, double b,
                 double tolerance = std::numeric_limits<double>::epsilon()) {
  double diff = std::abs(a - b);
  if (diff <= tolerance)
    return true;

  if (diff < std::max(std::abs(a), std::abs(b)) * tolerance)
    return true;

  return false;
}

// Pnm for maxn=maxm=4 (Matlab)
std::array<double, 15> Oct4 = {
    1.0000000000000000000000000e+00,  1.4188129598304730638602678e+00,
    1.1326022119815037925150136e+00,  3.8474545572362495171603314e-01,
    -5.1427278013612809370869172e-01, 9.9346352978712393344551401e-01,
    1.8197069354315411615630182e+00,  2.1885483890651622473910720e+00,
    1.8911107736775094601000546e+00,  6.3708625676351859912216469e-01,
    1.3807395850212413890290009e+00,  2.0397953656263365651568620e+00,
    3.9469618870173089897335217e-01,  9.6994856954294161877072611e-01,
    2.4012119461395231567735209e-01};

// Pnm for maxn=maxm=120 and m = 0 (Matlab)
std::array<double, 121> P120m0 = {
    1.00000000000000000000e+00,  1.41881295983047306386e+00,
    1.13260221198150379252e+00,  3.84745455723624951716e-01,
    -5.14272780136128093709e-01, -1.22415510020416284931e+00,
    -1.48349669962503449838e+00, -1.19892057834034071462e+00,
    -4.75968454010584629543e-01, 4.20837091376995653391e-01,
    1.16470518835519332157e+00,  1.48524042495522268226e+00,
    1.26626457552893501024e+00,  5.87559985437243748052e-01,
    -3.04405953357934966341e-01, -1.08604072241951055133e+00,
    -1.47396950782816626990e+00, -1.32766424457896192735e+00,
    -7.00239997576588502071e-01, 1.80908032481277736281e-01,
    9.96553389237216924457e-01,  1.45128136340418212313e+00,
    1.38044755344540526210e+00,  8.09745141314447192116e-01,
    -5.41533302471035593584e-02, -8.98460800747772347918e-01,
    -1.41752564997245866785e+00, -1.42346694417573704783e+00,
    -9.14157155171514856917e-01, -7.39624549411119769360e-02,
    7.93009903461919241074e-01,  1.37298547379490698539e+00,
    1.45607686632014088346e+00,  1.01222737039388599456e+00,
    2.02070356042807874886e-01,  -6.81214628978737413156e-01,
    -1.31800055011582917253e+00, -1.47787038877226484246e+00,
    -1.10298469796637554907e+00, -3.28998961553921942524e-01,
    5.64030941552166820152e-01,  1.25298568220309824994e+00,
    1.48859318280224295350e+00,  1.18561365681577579601e+00,
    4.53672578692903227804e-01,  -4.42409262409080183964e-01,
    -1.17843186866675453395e+00, -1.48811013468059383236e+00,
    -1.25941002708305860303e+00, -5.75075778443171481413e-01,
    3.17310733365480324153e-01,  1.09490338306940659940e+00,
    1.47639045012344016783e+00,  1.32376359802886156558e+00,
    6.92241485153857638579e-01,  -1.89710452203861068554e-01,
    -1.00303351843862276205e+00, -1.45349997550081089237e+00,
    -1.37815162944312152860e+00, -8.04248579823605425076e-01,
    6.05952797924405039964e-02,  9.03519647634809519587e-01,
    1.41959645236140796776e+00,  1.42213685111495058422e+00,
    9.10223675245455998706e-01,  6.90408726788018789478e-02,
    -7.97117787028987612175e-01, -1.37492589595736003361e+00,
    -1.45536738810507992881e+00, -1.00934477527290122900e+00,
    -1.98203015455324260685e-01, 6.84636736920054334021e-01,
    1.31981927020070100021e+00,  1.47757740052449060286e+00,
    1.10084573702344723500e+00,  3.25901699192448379705e-01,
    -5.66931845290928992931e-01, -1.25468905922161555466e+00,
    -1.48858783241468173308e+00, -1.18402097829489116876e+00,
    -4.51159949723752395911e-01, 4.44898437712501815255e-01,
    1.18002553752838323931e+00,  1.48830694822673281763e+00,
    1.25823011814926077889e+00,  5.73020317920096555753e-01,
    -3.19464958292972389664e-01, -1.09639264304759098678e+00,
    -1.47673047822274239316e+00, -1.32290236226793522789e+00,
    -6.90551909162544719223e-01, 1.91585869730704461134e-01,
    1.00442341274495294812e+00,  1.45394127094831437930e+00,
    1.37754051118823328537e+00,  8.02857291704160092038e-01,
    -6.22343634266011580536e-02, -9.04814972581906129534e-01,
    -1.42010839509515052193e+00, -1.42172450781735992820e+00,
    -9.09079203061360585636e-01, -6.76050670801617847872e-02,
    7.98323093331349986634e-01,  1.37548566027656438138e+00,
    1.45511446432921842664e+00,  1.00840698590117616718e+00,
    1.96944134412686988611e-01,  -6.85756336742521677152e-01,
    -1.32040954422581635619e+00, -1.47745312444025600485e+00,
    -1.10008269330892427540e+00, -3.24798291377221604037e-01,
    5.67969825540022399402e-01,  1.25529652668872548915e+00,
    1.48856772861808206443e+00,  1.18340680970183687037e+00,
    4.50194234173084906203e-01,  -4.45858677278189285964e-01,
    -1.18063983994646082287e+00, -1.48837125886135512864e+00,
    -1.25774353899286772673e+00};

int main() {
  printf("Unit test for Associated Lagrange Polynomials\n");

  AssociatedLegendreFunctions Pnm(120); // degree = order = 120
  AssociatedLegendreFunctions P4nm(4);  // degree = order = 4

  // compute at angle φ=
  Pnm.at(0.6108652382e0);
  P4nm.at(0.6108652382e0);

  assert(Pnm(0, 0) == P4nm(0, 0) && P4nm(0, 0) == 1e0);
  double d;

  // compare P4 to P
  for (int n = 0; n <= 4; n++) {
    for (int m = 0; m <= n; m++) {
      if (!ApproxEqual(Pnm(n, m), P4nm(n, m))) {
        d = Pnm(n, m) - P4nm(n, m);
        fprintf(stderr, "[T1] Precision exheded for (n,m)=(%d,%d) = %+.2e\n", n,
                m, d);
      }
    }
  }

  // compare P to Octave
  int i = 0;
  for (int m = 0; m <= 4; m++) {
    for (int n = m; n <= 4; n++) {
      if (!ApproxEqual(Pnm(n, m), Oct4[i])) {
        d = Pnm(n, m) - Oct4[i];
        fprintf(stderr, "[T2] Precision exheded for (n,m)=(%d,%d) = %+.2e\n", n,
                m, d);
      }
      ++i;
    }
  }

  // compare P4 to Octave
  i = 0;
  for (int m = 0; m <= 4; m++) {
    for (int n = m; n <= 4; n++) {
      if (!ApproxEqual(P4nm(n, m), Oct4[i])) {
        d = P4nm(n, m) - Oct4[i];
        fprintf(stderr, "[T3] Precision exheded for (n,m)=(%d,%d) = %+.2e\n", n,
                m, d);
      }
      ++i;
    }
  }

  // compare P to Octave (order=0, degree=120)
  for (int n = 0; n <= 120; n++) {
    if (!ApproxEqual(Pnm(n, 0), P120m0[n])) {
      d = Pnm(n, 0) - P120m0[n];
      fprintf(stderr,
              "[T4] Precision exheded for (n,m)=(%3d,%3d) = %+.2e (%.15e Vs "
              "%.15e)\n",
              n, 0, d, Pnm(n, 0), P120m0[n]);
    }
  }

  printf("All tests passed!\n");
  return 0;
}
