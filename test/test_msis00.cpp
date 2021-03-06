#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "atmosphere.hpp"
#include "doctest.h"
#include <cassert>
#include <cstdio>

using dso::air_density_models::nrlmsise00::gtd7;
#define DPRECISION 1e-12
#define TPRECISION 1e-9
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

TEST_CASE("MSIS00 Comparisson of results based on C implementation.\nSee "
          "https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/"
          "nrlmsis00_c_version/nrlmsise-00_test.c") {
  // input values
  double apha[7];
  for (int i = 0; i < 7; i++)
    apha[i] = 100;

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
  for (int i = 0; i < 15; i++) {
    gtd7(switches, input[i].doy, input[i].sec, input[i].lat, input[i].lon,
         input[i].lst, input[i].f107, input[i].f107A, input[i].alt, input[i].ap,
         apha, out[i].d, out[i].t);
  }

  // evaluate 15 and 16
  switches[9] = -1;
  for (int i = 15; i < 17; i++) {
    gtd7(switches, input[i].doy, input[i].sec, input[i].lat, input[i].lon,
         input[i].lst, input[i].f107, input[i].f107A, input[i].alt, input[i].ap,
         apha, out[i].d, out[i].t);
  }

  SUBCASE("MSIS00 Comparisson of results based on C implementation #0") {
    REQUIRE(std::abs(out[0].d[0] - 666517.69049515202641487) < DPRECISION);
    REQUIRE(std::abs(out[0].d[1] - 113880555.97522167861461639) < DPRECISION);
    REQUIRE(std::abs(out[0].d[2] - 19982109.25573454424738884) < DPRECISION);
    REQUIRE(std::abs(out[0].d[3] - 402276.35857125109760091) < DPRECISION);
    REQUIRE(std::abs(out[0].d[4] - 3557.46499451588579177) < DPRECISION);
    REQUIRE(std::abs(out[0].d[5] - 0.00000000000000407) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[0].d[6] - 34753.12399717141670408) < DPRECISION);
    REQUIRE(std::abs(out[0].d[7] - 4095913.26829300168901682) < DPRECISION);
    REQUIRE(std::abs(out[0].d[8] - 26672.73209335868887138) < DPRECISION);
    REQUIRE(std::abs(out[0].t[0] - 1250.53994356079942918) < TPRECISION);
    REQUIRE(std::abs(out[0].t[1] - 1241.41613001912060099) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #1") {
    REQUIRE(std::abs(out[1].d[0] - 3407293.22316091321408749) < DPRECISION);
    REQUIRE(std::abs(out[1].d[1] - 158633336.95691680908203125) < DPRECISION);
    REQUIRE(std::abs(out[1].d[2] - 13911173.65461114980280399) < DPRECISION);
    REQUIRE(std::abs(out[1].d[3] - 326255.95095955464057624) < DPRECISION);
    REQUIRE(std::abs(out[1].d[4] - 1559.61815050122459070) < DPRECISION);
    REQUIRE(std::abs(out[1].d[5] - 0.00000000000000500) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[1].d[6] - 48542.08463340254093055) < DPRECISION);
    REQUIRE(std::abs(out[1].d[7] - 4380966.71289862506091595) < DPRECISION);
    REQUIRE(std::abs(out[1].d[8] - 6956.68195594226835965) < DPRECISION);
    REQUIRE(std::abs(out[1].t[0] - 1166.75438375720887052) < TPRECISION);
    REQUIRE(std::abs(out[1].t[1] - 1161.71045188704238171) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #2") {
    REQUIRE(std::abs(out[2].d[0] - 112376.72440379358886275) < DPRECISION);
    REQUIRE(std::abs(out[2].d[1] - 69341.30086760609992780) < DPRECISION);
    REQUIRE(std::abs(out[2].d[2] - 42.47105217477091088) < DPRECISION);
    REQUIRE(std::abs(out[2].d[3] - 0.13227501414749321) < DPRECISION);
    REQUIRE(std::abs(out[2].d[4] - 0.00002618848418232) < DPRECISION);
    REQUIRE(std::abs(out[2].d[5] - 0.00000000000000000) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[2].d[6] - 20167.49854321431485005) < DPRECISION);
    REQUIRE(std::abs(out[2].d[7] - 5741.25593414717968699) < DPRECISION);
    REQUIRE(std::abs(out[2].d[8] - 23743.94151989594320185) < DPRECISION);
    REQUIRE(std::abs(out[2].t[0] - 1239.89211171666534028) < TPRECISION);
    REQUIRE(std::abs(out[2].t[1] - 1239.89064013305892331) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #3") {
    REQUIRE(std::abs(out[3].d[0] - 54115543.79936674237251282) < DPRECISION);
    REQUIRE(std::abs(out[3].d[1] - 191889344393.93090820312500000) <
            DPRECISION);
    REQUIRE(std::abs(out[3].d[2] - 6115825598224.63378906250000000) <
            DPRECISION);
    REQUIRE(std::abs(out[3].d[3] - 1225201051740.12426757812500000) <
            DPRECISION);
    REQUIRE(std::abs(out[3].d[4] - 60232119730.84870910644531250) < DPRECISION);
    REQUIRE(std::abs(out[3].d[5] - 0.00000000035844263) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[3].d[6] - 10598796.97740540280938148) < DPRECISION);
    REQUIRE(std::abs(out[3].d[7] - 261573.66937051419517957) < DPRECISION);
    REQUIRE(std::abs(out[3].d[8] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[3].t[0] - 1027.31846489999998084) < TPRECISION);
    REQUIRE(std::abs(out[3].t[1] - 206.88777640360549981) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #4") {
    REQUIRE(std::abs(out[4].d[0] - 1851122.48619252769276500) < DPRECISION);
    REQUIRE(std::abs(out[4].d[1] - 147655483.79274618625640869) < DPRECISION);
    REQUIRE(std::abs(out[4].d[2] - 15793562.28264496102929115) < DPRECISION);
    REQUIRE(std::abs(out[4].d[3] - 263379.49773123138584197) < DPRECISION);
    REQUIRE(std::abs(out[4].d[4] - 1588.78139838393008176) < DPRECISION);
    REQUIRE(std::abs(out[4].d[5] - 0.00000000000000481) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[4].d[6] - 58161.66780787473544478) < DPRECISION);
    REQUIRE(std::abs(out[4].d[7] - 5478984.47906878869980574) < DPRECISION);
    REQUIRE(std::abs(out[4].d[8] - 1264.44594176100849836) < DPRECISION);
    REQUIRE(std::abs(out[4].t[0] - 1212.39615212120929755) < TPRECISION);
    REQUIRE(std::abs(out[4].t[1] - 1208.13542521239173766) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #5") {
    REQUIRE(std::abs(out[5].d[0] - 867309.52339061570819467) < DPRECISION);
    REQUIRE(std::abs(out[5].d[1] - 127886176.80141274631023407) < DPRECISION);
    REQUIRE(std::abs(out[5].d[2] - 18225766.27171700075268745) < DPRECISION);
    REQUIRE(std::abs(out[5].d[3] - 292221.41906182467937469) < DPRECISION);
    REQUIRE(std::abs(out[5].d[4] - 2402.96243642370063753) < DPRECISION);
    REQUIRE(std::abs(out[5].d[5] - 0.00000000000000436) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[5].d[6] - 36863.89243750549940160) < DPRECISION);
    REQUIRE(std::abs(out[5].d[7] - 3897275.50372696388512850) < DPRECISION);
    REQUIRE(std::abs(out[5].d[8] - 26672.73209335868887138) < DPRECISION);
    REQUIRE(std::abs(out[5].t[0] - 1220.14641791503208879) < TPRECISION);
    REQUIRE(std::abs(out[5].t[1] - 1212.71208321180620260) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #6") {
    REQUIRE(std::abs(out[6].d[0] - 577625.12160232441965491) < DPRECISION);
    REQUIRE(std::abs(out[6].d[1] - 69791386.93660195171833038) < DPRECISION);
    REQUIRE(std::abs(out[6].d[2] - 12368135.59821702726185322) < DPRECISION);
    REQUIRE(std::abs(out[6].d[3] - 249286.77154291022452526) < DPRECISION);
    REQUIRE(std::abs(out[6].d[4] - 1405.73867417784299505) < DPRECISION);
    REQUIRE(std::abs(out[6].d[5] - 0.00000000000000247) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[6].d[6] - 52919.85567066640942357) < DPRECISION);
    REQUIRE(std::abs(out[6].d[7] - 1069814.10936656617559493) < DPRECISION);
    REQUIRE(std::abs(out[6].d[8] - 26672.73209335867795744) < DPRECISION);
    REQUIRE(std::abs(out[6].t[0] - 1116.38537604315160934) < TPRECISION);
    REQUIRE(std::abs(out[6].t[1] - 1112.99856821731100354) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #7") {
    REQUIRE(std::abs(out[7].d[0] - 374030.41055076656630263) < DPRECISION);
    REQUIRE(std::abs(out[7].d[1] - 47827201.23611342161893845) < DPRECISION);
    REQUIRE(std::abs(out[7].d[2] - 5240380.03332420438528061) < DPRECISION);
    REQUIRE(std::abs(out[7].d[3] - 175987.46403906072373502) < DPRECISION);
    REQUIRE(std::abs(out[7].d[4] - 550.16487795699640628) < DPRECISION);
    REQUIRE(std::abs(out[7].d[5] - 0.00000000000000157) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[7].d[6] - 88967.75722935037629213) < DPRECISION);
    REQUIRE(std::abs(out[7].d[7] - 1979740.83623295510187745) < DPRECISION);
    REQUIRE(std::abs(out[7].d[8] - 9121.81487599149295420) < DPRECISION);
    REQUIRE(std::abs(out[7].t[0] - 1031.24744071455893391) < TPRECISION);
    REQUIRE(std::abs(out[7].t[1] - 1024.84849221300896716) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #8") {
    REQUIRE(std::abs(out[8].d[0] - 674833.87666236236691475) < DPRECISION);
    REQUIRE(std::abs(out[8].d[1] - 124531526.04437313973903656) < DPRECISION);
    REQUIRE(std::abs(out[8].d[2] - 23690095.41052985191345215) < DPRECISION);
    REQUIRE(std::abs(out[8].d[3] - 491158.31547498231520876) < DPRECISION);
    REQUIRE(std::abs(out[8].d[4] - 4578.78109905442033778) < DPRECISION);
    REQUIRE(std::abs(out[8].d[5] - 0.00000000000000456) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[8].d[6] - 32445.94775161093275528) < DPRECISION);
    REQUIRE(std::abs(out[8].d[7] - 5370833.08708603773266077) < DPRECISION);
    REQUIRE(std::abs(out[8].d[8] - 26672.73209335868887138) < DPRECISION);
    REQUIRE(std::abs(out[8].t[0] - 1306.05204202729214558) < TPRECISION);
    REQUIRE(std::abs(out[8].t[1] - 1293.37404038953400232) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #9") {
    REQUIRE(std::abs(out[9].d[0] - 552860.08416451886296272) < DPRECISION);
    REQUIRE(std::abs(out[9].d[1] - 119804132.40413580834865570) < DPRECISION);
    REQUIRE(std::abs(out[9].d[2] - 34957977.64558214694261551) < DPRECISION);
    REQUIRE(std::abs(out[9].d[3] - 933961.83550281578209251) < DPRECISION);
    REQUIRE(std::abs(out[9].d[4] - 10962.54765493430204515) < DPRECISION);
    REQUIRE(std::abs(out[9].d[5] - 0.00000000000000497) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[9].d[6] - 26864.27856259808686445) < DPRECISION);
    REQUIRE(std::abs(out[9].d[7] - 4889974.23297140188515186) < DPRECISION);
    REQUIRE(std::abs(out[9].d[8] - 28054.44837125663616462) < DPRECISION);
    REQUIRE(std::abs(out[9].t[0] - 1361.86802078492337387) < TPRECISION);
    REQUIRE(std::abs(out[9].t[1] - 1347.38918372970169912) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #10") {
    REQUIRE(std::abs(out[10].d[0] - 137548758418628.65625000000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[10].d[1] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[10].d[2] - 20496870442907566080.00000000000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[10].d[3] - 5498695433718812672.00000000000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[10].d[4] - 245173315802838848.00000000000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[10].d[5] - 0.00126106566111855) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[10].d[6] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[10].d[7] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[10].d[8] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[10].t[0] - 1027.31846489999998084) < TPRECISION);
    REQUIRE(std::abs(out[10].t[1] - 281.46475766321560741) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #11") {
    REQUIRE(std::abs(out[11].d[0] - 44274425876770.92968750000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[11].d[1] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[11].d[2] - 6597567157737311232.00000000000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[11].d[3] - 1769929341406188544.00000000000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[11].d[4] - 78916799557274848.00000000000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[11].d[5] - 0.00040591393757992) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[11].d[6] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[11].d[7] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[11].d[8] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[11].t[0] - 1027.31846489999998084) < TPRECISION);
    REQUIRE(std::abs(out[11].t[1] - 227.41798082726180041) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #12") {
    REQUIRE(std::abs(out[12].d[0] - 2127828756207.18823242187500000) <
            DPRECISION);
    REQUIRE(std::abs(out[12].d[1] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[12].d[2] - 317079055035404288.00000000000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[12].d[3] - 85062798094347904.00000000000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[12].d[4] - 3792741116805988.50000000000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[12].d[5] - 0.00001950822245176) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[12].d[6] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[12].d[7] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[12].d[8] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[12].t[0] - 1027.31846489999998084) < TPRECISION);
    REQUIRE(std::abs(out[12].t[1] - 237.43891458772688452) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #13") {
    REQUIRE(std::abs(out[13].d[0] - 141218354559.28512573242187500) <
            DPRECISION);
    REQUIRE(std::abs(out[13].d[1] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[13].d[2] - 21043696437831580.00000000000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[13].d[3] - 5645392443377078.00000000000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[13].d[4] - 251714174941122.43750000000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[13].d[5] - 0.00000129470901593) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[13].d[6] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[13].d[7] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[13].d[8] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[13].t[0] - 1027.31846489999998084) < TPRECISION);
    REQUIRE(std::abs(out[13].t[1] - 279.55511295412799200) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #14") {
    REQUIRE(std::abs(out[14].d[0] - 12548844002.72669792175292969) <
            DPRECISION);
    REQUIRE(std::abs(out[14].d[1] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[14].d[2] - 1874532829219022.75000000000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[14].d[3] - 492305098078476.37500000000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[14].d[4] - 22396854138563.83984375000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[14].d[5] - 0.00000011476676715) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[14].d[6] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[14].d[7] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[14].d[8] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[14].t[0] - 1027.31846489999998084) < TPRECISION);
    REQUIRE(std::abs(out[14].t[1] - 219.07323136419572052) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #15") {
    REQUIRE(std::abs(out[15].d[0] - 519647.74029728799359873) < DPRECISION);
    REQUIRE(std::abs(out[15].d[1] - 127449407.29604625701904297) < DPRECISION);
    REQUIRE(std::abs(out[15].d[2] - 48504498.69853357225656509) < DPRECISION);
    REQUIRE(std::abs(out[15].d[3] - 1720837.98257490037940443) < DPRECISION);
    REQUIRE(std::abs(out[15].d[4] - 23544.86590544426144334) < DPRECISION);
    REQUIRE(std::abs(out[15].d[5] - 0.00000000000000588) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[15].d[6] - 25000.78391080929577583) < DPRECISION);
    REQUIRE(std::abs(out[15].d[7] - 6279209.82501879986375570) < DPRECISION);
    REQUIRE(std::abs(out[15].d[8] - 26672.73209335867795744) < DPRECISION);
    REQUIRE(std::abs(out[15].t[0] - 1426.41166228242468605) < TPRECISION);
    REQUIRE(std::abs(out[15].t[1] - 1408.60779555326394075) < TPRECISION);
  }
  SUBCASE("MSIS00 Comparisson of results based on C implementation #16") {
    REQUIRE(std::abs(out[16].d[0] - 42608597.48794120550155640) < DPRECISION);
    REQUIRE(std::abs(out[16].d[1] - 124134201554.87429809570312500) <
            DPRECISION);
    REQUIRE(std::abs(out[16].d[2] - 4929561542488.14257812500000000) <
            DPRECISION);
    REQUIRE(std::abs(out[16].d[3] - 1048406749092.83203125000000000) <
            DPRECISION);
    REQUIRE(std::abs(out[16].d[4] - 49934650830.55505371093750000) <
            DPRECISION);
    REQUIRE(std::abs(out[16].d[5] - 0.00000000029143036) <
            TOTALD_PRECISION); // total mass density
    REQUIRE(std::abs(out[16].d[6] - 8831228.59257159382104874) < DPRECISION);
    REQUIRE(std::abs(out[16].d[7] - 225251.55086261546239257) < DPRECISION);
    REQUIRE(std::abs(out[16].d[8] - 0.00000000000000000) < DPRECISION);
    REQUIRE(std::abs(out[16].t[0] - 1027.31846489999998084) < TPRECISION);
    REQUIRE(std::abs(out[16].t[1] - 193.40710625766814701) < TPRECISION);
  }
}
