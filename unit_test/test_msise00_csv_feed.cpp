#include "atmosphere.hpp"
#include "datetime/datetime_write.hpp" // strftime_ymd_hmfs
#include <cmath>
#include <cstdio>

constexpr const int num_tests = 10;

int differ(dso::nrlmsise00::detail::InParamsCore &in1,
           dso::nrlmsise00::detail::InParamsCore &in2) {
  int error = 0;
  if (std::abs(in1.year - in2.year) > 1e-6) {
    printf("\t> note failed on year: %d - %d\n", in1.year, in2.year);
    ++error;
  }
  if (std::abs(in1.doy - in2.doy) > 1e-6) {
    printf("\t> note failed on doy: %d - %d\n", in1.doy, in2.doy);
    ++error;
  }
  if (std::abs(in1.sec - in2.sec) > 1e-6) {
    printf("\t> note failed on sec: %.6f - %.6f\n", in1.sec, in2.sec);
    ++error;
  }
  if (std::abs(in1.f107A - in2.f107A) > 1e-6) {
    printf("\t> note failed on f107A: %.6f - %.6f\n", in1.f107A, in2.f107A);
    ++error;
  }
  if (std::abs(in1.f107 - in2.f107) > 1e-6) {
    printf("\t> note failed on f107: %.6f - %.6f\n", in1.f107, in2.f107);
    ++error;
  }
  if (std::abs(in1.ap - in2.ap) > 1e-6) {
    printf("\t> note failed on ap: %.6f - %.6f\n", in1.ap, in2.ap);
    ++error;
  }
  for (int i = 0; i < 7; i++) {
    if (std::abs(in1.aparr.a[i] - in2.aparr.a[i]) > 1e-6) {
      printf("\t> note failed on aparr.a[i]: %.6f - %.6f\n", in1.aparr.a[i],
             in2.aparr.a[i]);
      ++error;
    }
  }
  return error;
}

void print_params(dso::nrlmsise00::detail::InParamsCore &in,
                  const char *title) {
  printf("\t%s // %4d%2d:%.9f f107A: %.3f f107: %.3f Ap: %.3f [", title,
         in.year, in.doy, in.sec, in.f107A, in.f107, in.ap);
  for (int i = 0; i < 7; i++)
    printf("%.3f ", in.aparr.a[i]);
  printf("]\n");
  return;
}

const double sec_per_step = 60 * 60; // aka, 1 hour
const double sec_initial = 0e0;

dso::nrlmsise00::InParams<dso::nrlmsise00::detail::FluxDataFeedType::NONE>
    Results[num_tests];
/*
#
--------------------------------------------------------------------------------------------------------------------------------
#                              SPACE WEATHER DATA
#
--------------------------------------------------------------------------------------------------------------------------------
#
# See http://celestrak.com/SpaceData/SpaceWx-format.asp for format details.
#
# FORMAT(I4,I3,I3,I5,I3,8I3,I4,8I4,I4,F4.1,I2,I4,F6.1,I2,5F6.1)
#
# --------------------------------------------------------------------------------------------------------------------------------
#                                                                                             Adj     Adj   Adj   Obs   Obs   Obs 
# yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Q Ctr81 Lst81 F10.7 Ctr81 Lst81
# --------------------------------------------------------------------------------------------------------------------------------
#
2022 07 22 2577 10 37 23 27 23 13 23 17 27 190  22   9  12   9   5   9   6  12  10 0.6 3  99 118.4 0 120.7 132.3 114.7 117.2 128.6
2022 07 23 2577 11 43 57 37 20 20  7 13 27 224  32  67  22   7   7   3   5  12  19 1.0 5  97 114.1 0 120.6 132.3 110.5 117.0 128.6
2022 07 24 2577 12 27 23 20 13 13 23 23 17 159  12   9   7   5   5   9   9   6   7 0.4 2  79 110.5 0 120.3 132.0 107.1 116.8 128.3
2022 07 25 2577 13 27 17  3  7 10 10  7 20 101  12   6   2   3   4   4   3   7   5 0.2 1  86 105.5 0 119.8 131.8 102.3 116.4 128.1
2022 07 26 2577 14 33 20  7 17 17 10  7 20 131  18   7   3   6   6   4   3   7   6 0.3 1  91 101.9 0 119.5 131.6  98.8 116.0 127.8
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
    Results[i].params_.f107A = 116.4e0;
    Results[i].params_.f107 = 107.1e0;
    Results[i].params_.ap = 5e0;
    Results[i].params_.aparr = dso::nrlmsise00::detail::ApArray(
        {5e0, 12e0, 6e0, 9e0, 9e0,
         (double)(3 + 5 + 12 + 12 + 9 + 7 + 5 + 5) / 8,
         (double)(9 + 6 + 12 + 32 + 67 + 22 + 7 + 7) / 8e0});
  }
  // we will be adding 4 hours here, so it belongs to the next 3-hour interval
  for (int i = 3; i < 5; i++) {
    Results[i].params_.f107A = 116.4e0;
    Results[i].params_.f107 = 107.1e0;
    Results[i].params_.ap = 5e0;
    Results[i].params_.aparr = dso::nrlmsise00::detail::ApArray(
        {5e0, 6e0, 12e0, 6e0, 9e0,
         (double)(5 + 12 + 12 + 9 + 7 + 5 + 5 + 9) / 8,
         (double)(6 + 12 + 32 + 67 + 22 + 7 + 7 + 3) / 8e0});
  }
  // Done. let's try the example

  // 2022-07-25, doy=206
  dso::datetime<dso::seconds> t(dso::year(2022), dso::month(7),
                                dso::day_of_month(25), dso::seconds(0));

  dso::nrlmsise00::InParams<
      dso::nrlmsise00::detail::FluxDataFeedType::ST_CSV_SW>
      dataHunter(argv[1], t.mjd(), 0e0);

  int test_num = 0;
  int error = 0;
  char buf[128];

  // get results starting from 00 hours and adding one hour, a total of 5
  // times
  for (int step = 0; step < 5; step++) {
    dataHunter.update_params(t.mjd().as_underlying_type(), step * sec_per_step);
    const dso::datetime<dso::milliseconds> dt(
        t.mjd(), dso::milliseconds((long)(step * sec_per_step)));
    printf("> Test Case for %s\n", dso::strftime_ymd_hmfs(dt, buf));
    print_params(dataHunter.params_, "Test     ");
    print_params(Results[step].params_, "Reference");
    if (differ(dataHunter.params_, Results[step].params_)) {
      fprintf(stderr, "Failed! results differ for test case %d\n", step);
      ++error;
    }
    ++test_num;
  }

  // one day + 3 hrs forward ....
  // set the reference values
  Results[test_num].params_.year = 2022;
  Results[test_num].params_.doy = 206 + 1;
  Results[test_num].params_.sec = 3 * sec_per_step;
  Results[test_num].params_.f107A = 116.0;
  Results[test_num].params_.f107 = 102.3e0;
  Results[test_num].params_.ap = 6e0;
  Results[test_num].params_.aparr = dso::nrlmsise00::detail::ApArray(
      {6e0, 7e0, 18e0, 7e0, 3e0, (double)(9 + 6 + 12 + 6 + 2 + 3 + 4 + 4) / 8,
       (double)(5 + 12 + 12 + 9 + 7 + 5 + 5 + 9) / 8e0});
  // update the data feed ...
  dataHunter.update_params(t.mjd().as_underlying_type(),
                           (24 + 3) * sec_per_step);
  dso::datetime<dso::milliseconds> dt(
      t.mjd(), dso::milliseconds((long)((24 + 3) * sec_per_step)));
  printf("> Test Case for %s\n", dso::strftime_ymd_hmfs(dt, buf));
  print_params(dataHunter.params_, "Test     ");
  print_params(Results[test_num].params_, "Reference");
  if (differ(dataHunter.params_, Results[test_num].params_)) {
    fprintf(stderr, "Failed! results differ for test case %d\n", test_num);
    ++error;
  }
  ++test_num;

  // let's go back in time, request (again) for starting date, 1 second in
  Results[test_num].params_.year = Results[0].params_.year;
  Results[test_num].params_.doy = Results[0].params_.doy;
  Results[test_num].params_.sec = 1e0;
  Results[test_num].params_.f107A = Results[0].params_.f107A;
  Results[test_num].params_.f107 = Results[0].params_.f107;
  Results[test_num].params_.ap = Results[0].params_.ap;
  Results[test_num].params_.aparr = Results[0].params_.aparr;
  // update the data feed ...
  dataHunter.update_params(t.mjd().as_underlying_type(), 1e0);
  dt = dso::datetime<dso::milliseconds>(t.mjd(), dso::milliseconds(long(1)));
  printf("> Test Case for %s\n", dso::strftime_ymd_hmfs(dt, buf));
  print_params(dataHunter.params_, "Test     ");
  print_params(Results[test_num].params_, "Reference");
  if (differ(dataHunter.params_, Results[test_num].params_)) {
    fprintf(stderr, "Failed! results differ for test case %d\n", test_num);
    ++error;
  }
  ++test_num;
  
  // let's go just a little under one day forward, e.g. 86395 secs
  Results[test_num].params_.year = Results[0].params_.year;
  Results[test_num].params_.doy = Results[0].params_.doy;
  Results[test_num].params_.sec = 86395e0;
  Results[test_num].params_.f107A = Results[0].params_.f107A;
  Results[test_num].params_.f107 = Results[0].params_.f107;
  Results[test_num].params_.ap = Results[0].params_.ap;
  Results[test_num].params_.aparr = dso::nrlmsise00::detail::ApArray(
      {5e0, 7e0, 3e0, 4e0, 4e0, (double)(5+9+9+6+12+6+2+3) / 8,
       (double)(7+3+5+12+12+9+7+5) / 8e0});
  // update the data feed ...
  dataHunter.update_params(t.mjd().as_underlying_type(), 86395e0);
  dt = dso::datetime<dso::milliseconds>(t.mjd(), dso::milliseconds(long(86395e0 *1e3)));
  printf("> Test Case for %s\n", dso::strftime_ymd_hmfs(dt, buf));
  print_params(dataHunter.params_, "Test     ");
  print_params(Results[test_num].params_, "Reference");
  if (differ(dataHunter.params_, Results[test_num].params_)) {
    fprintf(stderr, "Failed! results differ for test case %d\n", test_num);
    ++error;
  }


  return error;
}
