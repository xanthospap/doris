#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "atmosphere.hpp"
#include <cstdio>

using dso::air_density_models::nrlmsise00::gtd7;
#define DPRECISION 1e-6
#define TPRECISION 1e-3
#define TOTALD_PRECISION 1e-15

struct auxin {
    double doy;
    int year;
    double sec;
    double alt;
    double lat;
    double lon;
    double lst;
    double f107;
    double f107A;
    double ap;
};

struct auxout {
    double d[10];
    double t[2];
};

TEST_CASE("MSIS00 Comparisson of results based on C implementation.\nSee https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/nrlmsis00_c_version/nrlmsise-00_test.c") {
  // input values
  double apha[7];
  for (int i = 0; i < 7; i++) apha[i] = 100;
  
    int switches[24];
  switches[0] = 0;
  for (int i = 1; i < 24; i++)
    switches[i] = 1;
  
  auxin input[17];
  auxout out[17];
  
  for (int i = 0; i < 17; i++) {
    input[i].doy = 172;
    input[i].year = 0; // without effect
    input[i].sec = 29000;
    input[i].alt = 400;
    input[i].lat = 60;
    input[i].lon = -70;
    input[i].lst = 16;
    input[i].f107A = 150;
    input[i].f107 = 150;
    input[i].ap = 4; // magnetic_index
  }

  input[1].doy = 81;
  input[2].sec = 75000;
  input[2].alt = 1000;
  input[3].alt = 100;
  input[10].alt = 0;
  input[11].alt = 10;
  input[12].alt = 30;
  input[13].alt = 50;
  input[14].alt = 70;
  input[16].alt = 100;
  input[4].lat = 0;
  input[5].lon = 0;
  input[6].lst = 4;
  input[7].f107A = 70;
  input[8].f107 = 180;
  input[9].ap = 40;
  // input[15].ap_a = &aph;
  // input[16].ap_a = &aph;

  // evaluate 0 to 14
  for (int i = 0; i < 15; i++)
    gtd7(switches, input[i].doy, input[i].sec,
         input[i].lat, input[i].lon, input[i].lst, input[i].f107, input[i].f107A,
         input[i].alt, input[i].ap, apha, out[i].d, out[i].t);

  // evaluate 15 and 16
  switches[9] = -1;
  for (int i = 15; i < 17; i++)
    gtd7(switches, input[i].doy, input[i].sec,
         input[i].lat, input[i].lon, input[i].lst, input[i].f107, input[i].f107A,
         input[i].alt, input[i].ap, apha, out[i].d, out[i].t);

SUBCASE("MSIS00 Comparisson of results based on C implementation #0") {
CHECK(std::abs(out[0].d[0]-6.665176904951520e+05)<DPRECISION);
CHECK(std::abs(out[0].d[1]-1.138805559752217e+08)<DPRECISION);
CHECK(std::abs(out[0].d[2]-1.998210925573454e+07)<DPRECISION);
CHECK(std::abs(out[0].d[3]-4.022763585712511e+05)<DPRECISION);
CHECK(std::abs(out[0].d[4]-3.557464994515886e+03)<DPRECISION);
CHECK(std::abs(out[0].d[5]-4.07471353275722235e-15)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[0].d[6]-3.475312399717142e+04)<DPRECISION);
CHECK(std::abs(out[0].d[7]-4.095913268293002e+06)<DPRECISION);
CHECK(std::abs(out[0].d[8]-2.667273209335869e+04)<DPRECISION);
CHECK(std::abs(out[0].t[0]-1.250539943560799e+03)<TPRECISION);
CHECK(std::abs(out[0].t[1]-1.241416130019121e+03)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #1") {
CHECK(std::abs(out[1].d[0]-3.407293223160913e+06)<DPRECISION);
CHECK(std::abs(out[1].d[1]-1.586333369569168e+08)<DPRECISION);
CHECK(std::abs(out[1].d[2]-1.391117365461115e+07)<DPRECISION);
CHECK(std::abs(out[1].d[3]-3.262559509595546e+05)<DPRECISION);
CHECK(std::abs(out[1].d[4]-1.559618150501225e+03)<DPRECISION);
CHECK(std::abs(out[1].d[5]-5.00184572907224415e-15)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[1].d[6]-4.854208463340254e+04)<DPRECISION);
CHECK(std::abs(out[1].d[7]-4.380966712898625e+06)<DPRECISION);
CHECK(std::abs(out[1].d[8]-6.956681955942268e+03)<DPRECISION);
CHECK(std::abs(out[1].t[0]-1.166754383757209e+03)<TPRECISION);
CHECK(std::abs(out[1].t[1]-1.161710451887042e+03)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #2") {
CHECK(std::abs(out[2].d[0]-1.123767244037936e+05)<DPRECISION);
CHECK(std::abs(out[2].d[1]-6.934130086760610e+04)<DPRECISION);
CHECK(std::abs(out[2].d[2]-4.247105217477091e+01)<DPRECISION);
CHECK(std::abs(out[2].d[3]-1.322750141474932e-01)<DPRECISION);
CHECK(std::abs(out[2].d[4]-2.618848418232185e-05)<DPRECISION);
CHECK(std::abs(out[2].d[5]-2.75677231926887374e-18)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[2].d[6]-2.016749854321431e+04)<DPRECISION);
CHECK(std::abs(out[2].d[7]-5.741255934147180e+03)<DPRECISION);
CHECK(std::abs(out[2].d[8]-2.374394151989594e+04)<DPRECISION);
CHECK(std::abs(out[2].t[0]-1.239892111716665e+03)<TPRECISION);
CHECK(std::abs(out[2].t[1]-1.239890640133059e+03)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #3") {
CHECK(std::abs(out[3].d[0]-5.411554379936674e+07)<DPRECISION);
CHECK(std::abs(out[3].d[1]-1.918893443939309e+11)<DPRECISION);
CHECK(std::abs(out[3].d[2]-6.115825598224634e+12)<DPRECISION);
CHECK(std::abs(out[3].d[3]-1.225201051740124e+12)<DPRECISION);
CHECK(std::abs(out[3].d[4]-6.023211973084871e+10)<DPRECISION);
CHECK(std::abs(out[3].d[5]-3.58442630411333433e-10)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[3].d[6]-1.059879697740540e+07)<DPRECISION);
CHECK(std::abs(out[3].d[7]-2.615736693705142e+05)<DPRECISION);
CHECK(std::abs(out[3].d[8]-2.819879355928334e-42)<DPRECISION);
CHECK(std::abs(out[3].t[0]-1.027318464900000e+03)<TPRECISION);
CHECK(std::abs(out[3].t[1]-2.068877764036055e+02)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #4") {
CHECK(std::abs(out[4].d[0]-1.851122486192528e+06)<DPRECISION);
CHECK(std::abs(out[4].d[1]-1.476554837927462e+08)<DPRECISION);
CHECK(std::abs(out[4].d[2]-1.579356228264496e+07)<DPRECISION);
CHECK(std::abs(out[4].d[3]-2.633794977312314e+05)<DPRECISION);
CHECK(std::abs(out[4].d[4]-1.588781398383930e+03)<DPRECISION);
CHECK(std::abs(out[4].d[5]-4.80963023940745184e-15)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[4].d[6]-5.816166780787474e+04)<DPRECISION);
CHECK(std::abs(out[4].d[7]-5.478984479068789e+06)<DPRECISION);
CHECK(std::abs(out[4].d[8]-1.264445941761008e+03)<DPRECISION);
CHECK(std::abs(out[4].t[0]-1.212396152121209e+03)<TPRECISION);
CHECK(std::abs(out[4].t[1]-1.208135425212392e+03)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #5") {
CHECK(std::abs(out[5].d[0]-8.673095233906157e+05)<DPRECISION);
CHECK(std::abs(out[5].d[1]-1.278861768014127e+08)<DPRECISION);
CHECK(std::abs(out[5].d[2]-1.822576627171700e+07)<DPRECISION);
CHECK(std::abs(out[5].d[3]-2.922214190618247e+05)<DPRECISION);
CHECK(std::abs(out[5].d[4]-2.402962436423701e+03)<DPRECISION);
CHECK(std::abs(out[5].d[5]-4.35586564264464624e-15)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[5].d[6]-3.686389243750550e+04)<DPRECISION);
CHECK(std::abs(out[5].d[7]-3.897275503726964e+06)<DPRECISION);
CHECK(std::abs(out[5].d[8]-2.667273209335869e+04)<DPRECISION);
CHECK(std::abs(out[5].t[0]-1.220146417915032e+03)<TPRECISION);
CHECK(std::abs(out[5].t[1]-1.212712083211806e+03)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #6") {
CHECK(std::abs(out[6].d[0]-5.776251216023244e+05)<DPRECISION);
CHECK(std::abs(out[6].d[1]-6.979138693660195e+07)<DPRECISION);
CHECK(std::abs(out[6].d[2]-1.236813559821703e+07)<DPRECISION);
CHECK(std::abs(out[6].d[3]-2.492867715429102e+05)<DPRECISION);
CHECK(std::abs(out[6].d[4]-1.405738674177843e+03)<DPRECISION);
CHECK(std::abs(out[6].d[5]-2.47065139166313155e-15)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[6].d[6]-5.291985567066641e+04)<DPRECISION);
CHECK(std::abs(out[6].d[7]-1.069814109366566e+06)<DPRECISION);
CHECK(std::abs(out[6].d[8]-2.667273209335868e+04)<DPRECISION);
CHECK(std::abs(out[6].t[0]-1.116385376043152e+03)<TPRECISION);
CHECK(std::abs(out[6].t[1]-1.112998568217311e+03)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #7") {
CHECK(std::abs(out[7].d[0]-3.740304105507666e+05)<DPRECISION);
CHECK(std::abs(out[7].d[1]-4.782720123611342e+07)<DPRECISION);
CHECK(std::abs(out[7].d[2]-5.240380033324204e+06)<DPRECISION);
CHECK(std::abs(out[7].d[3]-1.759874640390607e+05)<DPRECISION);
CHECK(std::abs(out[7].d[4]-5.501648779569964e+02)<DPRECISION);
CHECK(std::abs(out[7].d[5]-1.57188873925484437e-15)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[7].d[6]-8.896775722935038e+04)<DPRECISION);
CHECK(std::abs(out[7].d[7]-1.979740836232955e+06)<DPRECISION);
CHECK(std::abs(out[7].d[8]-9.121814875991493e+03)<DPRECISION);
CHECK(std::abs(out[7].t[0]-1.031247440714559e+03)<TPRECISION);
CHECK(std::abs(out[7].t[1]-1.024848492213009e+03)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #8") {
CHECK(std::abs(out[8].d[0]-6.748338766623624e+05)<DPRECISION);
CHECK(std::abs(out[8].d[1]-1.245315260443731e+08)<DPRECISION);
CHECK(std::abs(out[8].d[2]-2.369009541052985e+07)<DPRECISION);
CHECK(std::abs(out[8].d[3]-4.911583154749823e+05)<DPRECISION);
CHECK(std::abs(out[8].d[4]-4.578781099054420e+03)<DPRECISION);
CHECK(std::abs(out[8].d[5]-4.56442024536117137e-15)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[8].d[6]-3.244594775161093e+04)<DPRECISION);
CHECK(std::abs(out[8].d[7]-5.370833087086038e+06)<DPRECISION);
CHECK(std::abs(out[8].d[8]-2.667273209335869e+04)<DPRECISION);
CHECK(std::abs(out[8].t[0]-1.306052042027292e+03)<TPRECISION);
CHECK(std::abs(out[8].t[1]-1.293374040389534e+03)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #9") {
CHECK(std::abs(out[9].d[0]-5.528600841645189e+05)<DPRECISION);
CHECK(std::abs(out[9].d[1]-1.198041324041358e+08)<DPRECISION);
CHECK(std::abs(out[9].d[2]-3.495797764558215e+07)<DPRECISION);
CHECK(std::abs(out[9].d[3]-9.339618355028158e+05)<DPRECISION);
CHECK(std::abs(out[9].d[4]-1.096254765493430e+04)<DPRECISION);
CHECK(std::abs(out[9].d[5]-4.97454311032222523e-15)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[9].d[6]-2.686427856259809e+04)<DPRECISION);
CHECK(std::abs(out[9].d[7]-4.889974232971402e+06)<DPRECISION);
CHECK(std::abs(out[9].d[8]-2.805444837125664e+04)<DPRECISION);
CHECK(std::abs(out[9].t[0]-1.361868020784923e+03)<TPRECISION);
CHECK(std::abs(out[9].t[1]-1.347389183729702e+03)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #10") {
CHECK(std::abs(out[10].d[0]-1.375487584186287e+14)<DPRECISION);
CHECK(std::abs(out[10].d[1]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[10].d[2]-2.049687044290757e+19)<DPRECISION);
CHECK(std::abs(out[10].d[3]-5.498695433718813e+18)<DPRECISION);
CHECK(std::abs(out[10].d[4]-2.451733158028388e+17)<DPRECISION);
CHECK(std::abs(out[10].d[5]-1.26106566111855141e-03)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[10].d[6]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[10].d[7]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[10].d[8]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[10].t[0]-1.027318464900000e+03)<TPRECISION);
CHECK(std::abs(out[10].t[1]-2.814647576632156e+02)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #11") {
CHECK(std::abs(out[11].d[0]-4.427442587677093e+13)<DPRECISION);
CHECK(std::abs(out[11].d[1]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[11].d[2]-6.597567157737311e+18)<DPRECISION);
CHECK(std::abs(out[11].d[3]-1.769929341406189e+18)<DPRECISION);
CHECK(std::abs(out[11].d[4]-7.891679955727485e+16)<DPRECISION);
CHECK(std::abs(out[11].d[5]-4.05913937579917880e-04)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[11].d[6]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[11].d[7]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[11].d[8]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[11].t[0]-1.027318464900000e+03)<TPRECISION);
CHECK(std::abs(out[11].t[1]-2.274179808272618e+02)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #12") {
CHECK(std::abs(out[12].d[0]-2.127828756207188e+12)<DPRECISION);
CHECK(std::abs(out[12].d[1]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[12].d[2]-3.170790550354043e+17)<DPRECISION);
CHECK(std::abs(out[12].d[3]-8.506279809434790e+16)<DPRECISION);
CHECK(std::abs(out[12].d[4]-3.792741116805988e+15)<DPRECISION);
CHECK(std::abs(out[12].d[5]-1.95082224517562141e-05)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[12].d[6]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[12].d[7]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[12].d[8]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[12].t[0]-1.027318464900000e+03)<TPRECISION);
CHECK(std::abs(out[12].t[1]-2.374389145877269e+02)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #13") {
CHECK(std::abs(out[13].d[0]-1.412183545592851e+11)<DPRECISION);
CHECK(std::abs(out[13].d[1]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[13].d[2]-2.104369643783158e+16)<DPRECISION);
CHECK(std::abs(out[13].d[3]-5.645392443377078e+15)<DPRECISION);
CHECK(std::abs(out[13].d[4]-2.517141749411224e+14)<DPRECISION);
CHECK(std::abs(out[13].d[5]-1.29470901592856713e-06)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[13].d[6]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[13].d[7]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[13].d[8]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[13].t[0]-1.027318464900000e+03)<TPRECISION);
CHECK(std::abs(out[13].t[1]-2.795551129541280e+02)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #14") {
CHECK(std::abs(out[14].d[0]-1.254884400272670e+10)<DPRECISION);
CHECK(std::abs(out[14].d[1]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[14].d[2]-1.874532829219023e+15)<DPRECISION);
CHECK(std::abs(out[14].d[3]-4.923050980784764e+14)<DPRECISION);
CHECK(std::abs(out[14].d[4]-2.239685413856384e+13)<DPRECISION);
CHECK(std::abs(out[14].d[5]-1.14766767151153651e-07)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[14].d[6]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[14].d[7]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[14].d[8]-0.000000000000000e+00)<DPRECISION);
CHECK(std::abs(out[14].t[0]-1.027318464900000e+03)<TPRECISION);
CHECK(std::abs(out[14].t[1]-2.190732313641957e+02)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #15") {
CHECK(std::abs(out[15].d[0]-5.196477402972880e+05)<DPRECISION);
CHECK(std::abs(out[15].d[1]-1.274494072960463e+08)<DPRECISION);
CHECK(std::abs(out[15].d[2]-4.850449869853357e+07)<DPRECISION);
CHECK(std::abs(out[15].d[3]-1.720837982574900e+06)<DPRECISION);
CHECK(std::abs(out[15].d[4]-2.354486590544426e+04)<DPRECISION);
CHECK(std::abs(out[15].d[5]-5.88194044865163181e-15)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[15].d[6]-2.500078391080930e+04)<DPRECISION);
CHECK(std::abs(out[15].d[7]-6.279209825018800e+06)<DPRECISION);
CHECK(std::abs(out[15].d[8]-2.667273209335868e+04)<DPRECISION);
CHECK(std::abs(out[15].t[0]-1.426411662282425e+03)<TPRECISION);
CHECK(std::abs(out[15].t[1]-1.408607795553264e+03)<TPRECISION);
}SUBCASE("MSIS00 Comparisson of results based on C implementation #16") {
CHECK(std::abs(out[16].d[0]-4.260859748794121e+07)<DPRECISION);
CHECK(std::abs(out[16].d[1]-1.241342015548743e+11)<DPRECISION);
CHECK(std::abs(out[16].d[2]-4.929561542488143e+12)<DPRECISION);
CHECK(std::abs(out[16].d[3]-1.048406749092832e+12)<DPRECISION);
CHECK(std::abs(out[16].d[4]-4.993465083055505e+10)<DPRECISION);
CHECK(std::abs(out[16].d[5]-2.91430355030879247e-10)<TOTALD_PRECISION);// total mass density
CHECK(std::abs(out[16].d[6]-8.831228592571594e+06)<DPRECISION);
CHECK(std::abs(out[16].d[7]-2.252515508626155e+05)<DPRECISION);
CHECK(std::abs(out[16].d[8]-2.415245929648914e-42)<DPRECISION);
CHECK(std::abs(out[16].t[0]-1.027318464900000e+03)<TPRECISION);
CHECK(std::abs(out[16].t[1]-1.934071062576681e+02)<TPRECISION);
}
}