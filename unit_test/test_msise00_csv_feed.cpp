#include "atmosphere.hpp"
#include <cmath>
#include <cstdio>

constexpr const int num_tests = 6;

int differ(dso::nrlmsise00::detail::InParamsCore &in1,
           dso::nrlmsise00::detail::InParamsCore &in2) {
  int error = 0;
  if (std::abs(in1.year - in2.year) > 1e-6)
    ++error;
  if (std::abs(in1.doy - in2.doy) > 1e-6)
    ++error;
  if (std::abs(in1.sec - in2.sec) > 1e-6)
    ++error;
  if (std::abs(in1.f107A - in2.f107A) > 1e-6)
    ++error;
  if (std::abs(in1.f107 - in2.f107) > 1e-6)
    ++error;
  if (std::abs(in1.ap - in2.ap) > 1e-6)
    ++error;
  for (int i = 0; i < 7; i++)
    if (std::abs(in1.aparr.a[i] - in2.aparr.a[i]) > 1e-6)
      ++error;
  return error;
}

void print_params(dso::nrlmsise00::detail::InParamsCore &in) {
  printf("\t%4d%2d:%.9f f107A: %.3f f107: %.3f Ap: %.3f [", in.year, in.doy,
         in.sec, in.f107A, in.f107, in.ap);
  for (int i = 0; i < 7; i++)
    printf("%.3f ", in.aparr.a[i]);
  printf("]\n");
  return;
}

const double sec_per_step = 60 * 60; // aka, 1 hour
const double sec_initial = 0e0;

dso::nrlmsise00::InParams<
  dso::nrlmsise00::detail::FluxDataFeedType::NONE>
    Results[num_tests];
/*
# --------------------------------------------------------------------------------------------------------------------------------
#                              SPACE WEATHER DATA
# --------------------------------------------------------------------------------------------------------------------------------
#
# See http://celestrak.com/SpaceData/SpaceWx-format.asp for format details.
#
# FORMAT(I4,I3,I3,I5,I3,8I3,I4,8I4,I4,F4.1,I2,I4,F6.1,I2,5F6.1)
# --------------------------------------------------------------------------------------------------------------------------------
#                                                                                             Adj     Adj   Adj   Obs   Obs   Obs 
# yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Q Ctr81 Lst81 F10.7 Ctr81 Lst81
# --------------------------------------------------------------------------------------------------------------------------------
2022 07 22 2577 10 37 23 27 23 13 23 17 27 190  22   9  12   9   5   9   6  12  10 0.6 3 103 118.4 0 120.9 132.3 114.7 117.4 128.6
2022 07 23 2577 11 43 57 37 20 20  7 13 27 224  32  67  22   7   7   3   5  12  19 1.0 5 102 114.1 0 120.8 132.3 110.5 117.3 128.6
2022 07 24 2577 12 27 23 20 13 13 23 23 17 159  12   9   7   5   5   9   9   6   7 0.4 2  84 110.5 0 120.5 132.0 107.1 117.0 128.3
2022 07 25 2577 13 27 17  3  7 10 10  7 20 101  12   6   2   3   4   4   3   7   5 0.2 1  87 105.5 0 120.1 131.8 102.3 116.6 128.1
2022 07 26 2577 14 33 20  7 17 17 10  7 20 131  18   7   3   6   6   4   3   7   6 0.3 1  83 101.9 0 119.7 131.6  98.8 116.3 127.8
*/

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr,
            "ERROR. Must provide CelesTrak Space-Weather CSV data file\n");
    return 1;
  }

// Fill in results ...
// Reference results for date/time
for (int i = 0; i < 5; i++) {
  Results[i].params_.year = 2022;
  Results[i].params_.doy = 206;
  Results[i].params_.sec = sec_initial + i * sec_per_step;
}

// reference results for flux/indexes
for (int i = 0; i < 3; i++) {
  Results[i].params_.f107A = 116.6;
  Results[i].params_.f107 = 107.1e0;
  Results[i].params_.ap = 5e0;
  Results[i].params_.aparr = dso::nrlmsise00::detail::ApArray(
      {5e0, 12e0, 6e0, 9e0, 9e0, (double)(3 + 5 + 12 + 12 + 9 + 7 + 5 + 5) / 8,
       (double)(9 +6 +12 + 32 +67 +22 +7  +7) / 8e0});
}
// we will be adding 4 hours here, so it belongs to the next 3-hour interval
for (int i = 3; i < 5; i++) {
  Results[i].params_.f107A = 116.6;
  Results[i].params_.f107 = 107.1e0;
  Results[i].params_.ap = 5e0;
  Results[i].params_.aparr = dso::nrlmsise00::detail::ApArray(
      {5e0, 6e0, 12e0, 6e0, 9e0, (double)(5 + 12 + 12 + 9 + 7 + 5 + 5 + 9) / 8,
       (double)(6 +12 + 32 +67 +22 +7  +7 +3) / 8e0});
}
// sixth test case is one day plus 3 hrs ahead
Results[5].params_.year = 2022;
Results[5].params_.doy = 206;
Results[5].params_.sec = sec_initial + 24 * sec_per_step + 3 * sec_per_step;
Results[5].params_.f107A = 116.3;
Results[5].params_.f107 = 98.8e0;
Results[5].params_.ap = 6e0;
Results[5].params_.aparr = dso::nrlmsise00::detail::ApArray(
    {6e0,7e0,18e0,7e0,3e0, (double)(9+6+6+2+3+4+4) / 8,
     (double)(5+12+12+9+7+5+5+9) / 8e0});
// Done. let's try the example

  // 2022-07-25, doy=206
  dso::datetime<dso::seconds> t(dso::year(2022), dso::month(7),
                                dso::day_of_month(25), dso::seconds(0));

  dso::nrlmsise00::InParams<
      dso::nrlmsise00::detail::FluxDataFeedType::ST_CSV_SW>
      dataHunter(argv[1], t.mjd(), 0e0);

  // get reults starting from 00 hours and adding one hour
  for (int step = 0; step < 5; step++) {
    printf("--------- Hours in day: %.3f ------------------\n",
           step * sec_per_step / (60e0 * 60e0));
    dataHunter.update_params(t.mjd().as_underlying_type(), step * sec_per_step);
    print_params(dataHunter.params_);
    print_params(Results[step].params_);
    if (differ(dataHunter.params_, Results[step].params_)) {
      fprintf(stderr, "Failed! results differ for test case %d\n", step);
      return 1;
    }
  }
  
  // one day + 3 hrs forward ....
  dataHunter.update_params(t.mjd().as_underlying_type(), (24+3) * sec_per_step);
    print_params(dataHunter.params_);
    print_params(Results[5].params_);
    if (differ(dataHunter.params_, Results[5].params_)) {
      fprintf(stderr, "Failed! results differ for test case %d\n", 5);
      return 1;
    }


  return 0;
}
