#include "egravity.hpp"
#include "geodesy/units.hpp"
#include "tides.hpp"
#include <stdexcept>
#ifdef DEBUG
#include <cassert>
#endif

namespace {
/* γ coefficients, real and imaginary part, see IERS2010, Sec. 6.4 */
constexpr const double gamma_real = 0.6870e0;
constexpr const double gamma_imag = 0.0036e0;

/* load deformation coefficients extracted from COST-G,
 * see loadLoveNumbers_Gegout97.txt
 */
constexpr const double LoadLoveNumbers[360] = {
    -0.3054020195e+00, -0.1960294041e+00, -0.1336652689e+00, -0.1047066267e+00,
    -0.9033564429e-01, -0.8206984804e-01, -0.7655494644e-01, -0.7243844815e-01,
    -0.6913401466e-01, -0.6635869819e-01, -0.6395689877e-01, -0.6183296641e-01,
    -0.5992172201e-01, -0.5817772516e-01, -0.5656704205e-01, -0.5506474110e-01,
    -0.5365205058e-01, -0.5231479422e-01, -0.5104204279e-01, -0.4982562670e-01,
    -0.4865919131e-01, -0.4753764652e-01, -0.4645728556e-01, -0.4541483359e-01,
    -0.4440818528e-01, -0.4343487603e-01, -0.4249347823e-01, -0.4158232061e-01,
    -0.4070034557e-01, -0.3984645832e-01, -0.3901937759e-01, -0.3821827509e-01,
    -0.3744217649e-01, -0.3669034206e-01, -0.3596186474e-01, -0.3525594412e-01,
    -0.3457187576e-01, -0.3390882631e-01, -0.3326609909e-01, -0.3264299923e-01,
    -0.3203883174e-01, -0.3145293280e-01, -0.3088463896e-01, -0.3033334191e-01,
    -0.2979842108e-01, -0.2927929373e-01, -0.2877539004e-01, -0.2828615868e-01,
    -0.2781108192e-01, -0.2734963353e-01, -0.2690133637e-01, -0.2646570945e-01,
    -0.2604229418e-01, -0.2563066630e-01, -0.2523039452e-01, -0.2484107774e-01,
    -0.2446233127e-01, -0.2409377795e-01, -0.2373506405e-01, -0.2338584530e-01,
    -0.2304579309e-01, -0.2271459030e-01, -0.2239193633e-01, -0.2207753919e-01,
    -0.2177111937e-01, -0.2147240908e-01, -0.2118115170e-01, -0.2089710139e-01,
    -0.2062002059e-01, -0.2034968220e-01, -0.2008586807e-01, -0.1982836895e-01,
    -0.1957698319e-01, -0.1933151728e-01, -0.1909178544e-01, -0.1885760922e-01,
    -0.1862881706e-01, -0.1840524275e-01, -0.1818672781e-01, -0.1797312161e-01,
    -0.1776427277e-01, -0.1756004164e-01, -0.1736029175e-01, -0.1716489163e-01,
    -0.1697371498e-01, -0.1678663942e-01, -0.1660354744e-01, -0.1642432560e-01,
    -0.1624886478e-01, -0.1607705911e-01, -0.1590880688e-01, -0.1574400976e-01,
    -0.1558257306e-01, -0.1542440486e-01, -0.1526941673e-01, -0.1511752309e-01,
    -0.1496864146e-01, -0.1482269168e-01, -0.1467959656e-01, -0.1453928135e-01,
    -0.1440167387e-01, -0.1426670401e-01, -0.1413430411e-01, -0.1400440861e-01,
    -0.1387695416e-01, -0.1375187918e-01, -0.1362912416e-01, -0.1350863142e-01,
    -0.1339034516e-01, -0.1327421107e-01, -0.1316017669e-01, -0.1304819106e-01,
    -0.1293820488e-01, -0.1283017014e-01, -0.1272404038e-01, -0.1261977050e-01,
    -0.1251731678e-01, -0.1241663664e-01, -0.1231768886e-01, -0.1222043336e-01,
    -0.1212483127e-01, -0.1203084476e-01, -0.1193843710e-01, -0.1184757260e-01,
    -0.1175821983e-01, -0.1167033892e-01, -0.1158389997e-01, -0.1149887118e-01,
    -0.1141522155e-01, -0.1133292103e-01, -0.1125194019e-01, -0.1117225053e-01,
    -0.1109382428e-01, -0.1101663454e-01, -0.1094065488e-01, -0.1086585976e-01,
    -0.1079222424e-01, -0.1071972412e-01, -0.1064833570e-01, -0.1057803597e-01,
    -0.1050880250e-01, -0.1044061352e-01, -0.1037344765e-01, -0.1030728418e-01,
    -0.1024210288e-01, -0.1017788410e-01, -0.1011460854e-01, -0.1005225749e-01,
    -0.9990812687e-02, -0.9930256346e-02, -0.9870571027e-02, -0.9811739803e-02,
    -0.9753746115e-02, -0.9696573861e-02, -0.9640207245e-02, -0.9584630905e-02,
    -0.9529829829e-02, -0.9475790189e-02, -0.9422496073e-02, -0.9369934305e-02,
    -0.9318091242e-02, -0.9266953580e-02, -0.9216508278e-02, -0.9166742633e-02,
    -0.9117644217e-02, -0.9069200890e-02, -0.9021400853e-02, -0.8974232435e-02,
    -0.8927684326e-02, -0.8881745447e-02, -0.8836405029e-02, -0.8791653407e-02,
    -0.8747481667e-02, -0.8703874181e-02, -0.8660824198e-02, -0.8618321971e-02,
    -0.8576358051e-02, -0.8534923154e-02, -0.8494008256e-02, -0.8453604422e-02,
    -0.8413702992e-02, -0.8374295452e-02, -0.8335373508e-02, -0.8296928978e-02,
    -0.8258953895e-02, -0.8221440446e-02, -0.8184381010e-02, -0.8147768058e-02,
    -0.8111594272e-02, -0.8075852456e-02, -0.8040535598e-02, -0.8005636771e-02,
    -0.7971149225e-02, -0.7937066335e-02, -0.7903381636e-02, -0.7870088750e-02,
    -0.7837181442e-02, -0.7804653597e-02, -0.7772499254e-02, -0.7740712516e-02,
    -0.7709287634e-02, -0.7678218959e-02, -0.7647500968e-02, -0.7617128214e-02,
    -0.7587095379e-02, -0.7557397242e-02, -0.7528028665e-02, -0.7498984668e-02,
    -0.7470260263e-02, -0.7441850636e-02, -0.7413751027e-02, -0.7385956807e-02,
    -0.7358463369e-02, -0.7331266234e-02, -0.7304360991e-02, -0.7277743343e-02,
    -0.7251409005e-02, -0.7225353835e-02, -0.7199573711e-02, -0.7174064656e-02,
    -0.7148822685e-02, -0.7123843938e-02, -0.7099129566e-02, -0.7074666551e-02,
    -0.7050455306e-02, -0.7026492490e-02, -0.7002774540e-02, -0.6979298005e-02,
    -0.6956059431e-02, -0.6933055469e-02, -0.6910282814e-02, -0.6887738228e-02,
    -0.6865418504e-02, -0.6843320529e-02, -0.6821441187e-02, -0.6799777491e-02,
    -0.6778326425e-02, -0.6757085071e-02, -0.6736050549e-02, -0.6715220049e-02,
    -0.6694590762e-02, -0.6674159952e-02, -0.6653925199e-02, -0.6633883871e-02,
    -0.6614032563e-02, -0.6594369232e-02, -0.6574891356e-02, -0.6555596461e-02,
    -0.6536482218e-02, -0.6517546015e-02, -0.6498785614e-02, -0.6480198683e-02,
    -0.6461782835e-02, -0.6443535971e-02, -0.6425455846e-02, -0.6407540263e-02,
    -0.6389786969e-02, -0.6372194033e-02, -0.6354759310e-02, -0.6337480778e-02,
    -0.6320356275e-02, -0.6303384010e-02, -0.6286561980e-02, -0.6269888285e-02,
    -0.6253360871e-02, -0.6236978078e-02, -0.6220737994e-02, -0.6204638851e-02,
    -0.6188678733e-02, -0.6172856062e-02, -0.6157169047e-02, -0.6141616004e-02,
    -0.6126195168e-02, -0.6110905001e-02, -0.6095743854e-02, -0.6080714363e-02,
    -0.6065808661e-02, -0.6051025499e-02, -0.6036365150e-02, -0.6021826113e-02,
    -0.6007406839e-02, -0.5993105878e-02, -0.5978921847e-02, -0.5964853301e-02,
    -0.5950898825e-02, -0.5937057006e-02, -0.5923326579e-02, -0.5909706151e-02,
    -0.5896194416e-02, -0.5882789995e-02, -0.5869491718e-02, -0.5856298304e-02,
    -0.5843208448e-02, -0.5830220904e-02, -0.5817334546e-02, -0.5804548161e-02,
    -0.5791860577e-02, -0.5779270522e-02, -0.5766776986e-02, -0.5754378816e-02,
    -0.5742074899e-02, -0.5729864025e-02, -0.5717745257e-02, -0.5705717481e-02,
    -0.5693779662e-02, -0.5681930652e-02, -0.5670169557e-02, -0.5658495321e-02,
    -0.5646906981e-02, -0.5635403426e-02, -0.5623983856e-02, -0.5612647214e-02,
    -0.5601392619e-02, -0.5590219039e-02, -0.5579125640e-02, -0.5568111507e-02,
    -0.5557175755e-02, -0.5546317422e-02, -0.5535535712e-02, -0.5524829768e-02,
    -0.5514198739e-02, -0.5503641734e-02, -0.5493157962e-02, -0.5482746632e-02,
    -0.5472406939e-02, -0.5462138054e-02, -0.5451939170e-02, -0.5441809594e-02,
    -0.5431748532e-02, -0.5421755212e-02, -0.5411828857e-02, -0.5401975569e-02,
    -0.5392181719e-02, -0.5382453176e-02, -0.5372788210e-02, -0.5363186748e-02,
    -0.5353648059e-02, -0.5344171504e-02, -0.5334756310e-02, -0.5325401926e-02,
    -0.5316107680e-02, -0.5306872920e-02, -0.5297696936e-02, -0.5288579215e-02,
    -0.5279519432e-02, -0.5270516389e-02, -0.5261569543e-02, -0.5252678535e-02,
    -0.5243842709e-02, -0.5235061516e-02, -0.5226334273e-02, -0.5217660524e-02,
    -0.5209039665e-02, -0.5200471169e-02, -0.5191954431e-02, -0.5183488963e-02,
    -0.5175074221e-02, -0.5166709698e-02, -0.5158394787e-02, -0.5150129056e-02};
} // namespace

dso::OceanPoleTide::OceanPoleTide(int maxdegree, int maxorder, const char *fn,
                                  double GMearth, double Rearth)
    : max_degree(maxdegree), max_order(maxorder),
      Ar(maxdegree + 1, maxorder + 1), Ai(maxdegree + 1, maxorder + 1),
      Br(maxdegree + 1, maxorder + 1), Bi(maxdegree + 1, maxorder + 1),
      dCS(maxdegree, maxorder, GMearth,
          Rearth), /*Rn(new double[maxdegree + 1]),*/
      V(maxdegree + 3, maxdegree + 3), W(maxdegree + 3, maxdegree + 3) {
  if (dso::parse_desai_ocean_pole_tide_coeffs(fn, max_degree, max_order, Ar, Ai,
                                              Br, Bi)) {
    throw std::runtime_error("Failed constructing OceanPoleTide instance!\n");
  }
  /* constants */
  const double Omega = iers2010::OmegaEarth;
  const double ae = /*iesr2010::Re*/ Rearth;
  const double GM = /*iers2010::GMe*/ GMearth;
  const double ge = iers2010::ge;
  const double G = iers2010::G;
  const double rho = 1025e0; /* density of water in [kg/m^3] */

  const double fac = ((Omega * Omega * ae * ae) / GM) * (ae * ae) * (G / ge) *
                     (4e0 * dso::DPI * rho);
  Rn.reserve(max_degree + 1);
  for (int i = 0; i <= max_degree; i++) {
    // Rn[i] = fac * (1e0 + LoadLoveNumbers[i]) / (2 * i + 1e0);
    Rn.push_back(fac * (1e0 + LoadLoveNumbers[i]) / (2 * i + 1e0));
  }
}

// dso::OceanPoleTide::~OceanPoleTide() noexcept {
//   if (Rn) delete[] Rn;
// }

int dso::OceanPoleTide::operator()(const dso::TwoPartDate &ttmjd, double xp,
                                   double yp, int maxdegree,
                                   int maxorder) noexcept {
  /* secular pole in [mas] */
  double xs, ys;
  secular_pole(ttmjd, xs, ys);
  /* wobble in [as] */
  double m1, m2;
  wobble_components(xp, yp, xs, ys, m1, m2);
  /* wobble in radians */
  m1 = dso::sec2rad(m1);
  m2 = dso::sec2rad(m2);
  
  if (maxdegree < 0)
    maxdegree = this->max_degree;
  if (maxorder < 0)
    maxorder = this->max_order;

  /* perform computation column-wise, starting from bigger (n,m) values */
  for (int m = maxorder; m >= 0; m--) {
    for (int n = maxdegree; n >= m; n--) {
      dCS.C(n, m) = Ar(n, m) * (m1 * gamma_real + m2 * gamma_imag) +
                    Ai(n, m) * (m2 * gamma_real - m1 * gamma_imag);
      dCS.C(n, m) *= Rn[n];
      dCS.S(n, m) = Br(n, m) * (m1 * gamma_real + m2 * gamma_imag) +
                    Bi(n, m) * (m2 * gamma_real - m1 * gamma_imag);
      dCS.S(n, m) *= Rn[n];
    }
  }

  return 0;
}

int dso::OceanPoleTide::acceleration(const dso::TwoPartDate &mjdtt,
                                     double xp_sec, double yp_sec,
                                     const Eigen::Matrix<double, 3, 1> &rsat,
                                     Eigen::Matrix<double, 3, 1> &acc,
                                     Eigen::Matrix<double, 3, 3> &acc_gradient,
                                     int maxdegree, int maxorder) noexcept {
  if (maxdegree < 0)
    maxdegree = this->max_degree;
  if (maxorder < 0)
    maxorder = this->max_order;
  /* prepare Stokes coefficients */
  this->operator()(mjdtt, xp_sec, yp_sec, maxdegree, maxorder);
  /* compute acceleration at satellite position (ITRF, cartesian) */
  dso::gravity_acceleration(dCS, rsat, maxdegree, dCS.Re(), dCS.GM(), acc,
                            acc_gradient, &V, &W);

  return 0;
}