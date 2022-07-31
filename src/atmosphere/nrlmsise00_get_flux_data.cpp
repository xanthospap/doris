#include "nrlmsise00.hpp"
#include <stdexcept>
#ifdef DEBUG
#include <cstdio>
#endif

dso::nrlmsise00::detail::Nrlmsise00DataFeed::Nrlmsise00DataFeed(
    const char *fncsv, dso::nrlmsise00::detail::InParamsCore &in)
    : fncsv_(fncsv) {

  const dso::datetime<dso::milliseconds> t(dso::datetime<dso::milliseconds>(
      dso::year(in.year), dso::day_of_year(in.doy), dso::milliseconds(0)));
  if (dso::utils::celestrak::details::parse_csv_for_date(t.mjd(), fncsv,
                                                         flux_data_, 3, 0)) {
    throw std::runtime_error("Failed to initialize Nrlmsise00DataFeed\n");
  }

  mjd_ = t.mjd();

  if (init(in))
    throw std::runtime_error("Failed to initialize Nrlmsise00DataFeed\n");
}

int ihours(double secday) noexcept {
  constexpr const int sec_in_hour = 60 * 60;
  return static_cast<int>(secday) / sec_in_hour;
}

int dso::nrlmsise00::detail::Nrlmsise00DataFeed::init(
    dso::nrlmsise00::detail::InParamsCore &in) noexcept {
#ifdef DEBUG
      for (int i=0; i<4; i++) {
        printf("%d:", (int)flux_data_[i].mjd_.as_underlying_type());
        for (int p=0; p<8; p++)
          printf(" %2d", flux_data_[i].ApIndexes[p]);
        printf("\n");
      }
#endif

  // flux data for this (most recent) day
  const dso::utils::celestrak::details::CelestTrakSWFlux *cur = flux_data_ + 3;
#ifdef DEBUG
  const dso::utils::celestrak::details::CelestTrakSWFlux *today = flux_data_ + 3;
#endif

  // F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
  in.f107A = cur->f107ObsC81;

  // F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
  in.f107 = cur[-1].f107ObsC81;

  // AP - MAGNETIC INDEX(DAILY)
  in.ap = cur->ApDailyAverage;

  // ap array, contains:
  // (1) DAILY AP
  in.aparr.a[0] = cur->ApDailyAverage;

  // // (2) 3 HR AP INDEX FOR CURRENT TIME
  int index = ihours(in.sec) / 3;
#ifdef DEBUG
  assert(index >= 0 && index < 8);
#endif
  in.aparr.a[1] = cur->ApIndexes[index];

  // (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
  --index;
  if (index < 0) {
    --cur;
    index = 7;
  }
  in.aparr.a[2] = cur->ApIndexes[index];
  printf("\t 3-hr index =%3d\n", cur->ApIndexes[index]);

  // (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
  --index;
  if (index < 0) {
    --cur;
    index = 7;
  }
  in.aparr.a[3] = cur->ApIndexes[index];
  printf("\t 6-hr index =%3d\n", cur->ApIndexes[index]);

  // (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
  --index;
  if (index < 0) {
    --cur;
    index = 7;
  }
  in.aparr.a[4] = cur->ApIndexes[index];
  printf("\t 9-hr index =%3d\n", cur->ApIndexes[index]);

  // (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
  // TO CURRENT TIME
  int stop_index = index;
  auto stop_doy = cur;
#ifdef DEBUG
  printf("> Aps for -33 to -12 hours ...\n");
  printf("\tCurrent hour index = %d\n", ihours(in.sec) % 3);
#endif
  // now, 33 hours is 24 + 9
  int start_index = ihours(in.sec) / 3 - 3 + 1;
  const dso::utils::celestrak::details::CelestTrakSWFlux *start_doy = flux_data_ + 3 - 1;
  if (start_index < 0) {
    --start_doy;
    start_index += 7;
  }
#ifdef DEBUG
  printf("\tWill start at %ld days before and index %d\n", start_doy - today, start_index);
  printf("\tWill stop at %ld days before and index %d\n", stop_doy - today, stop_index);
  assert(start_index >= 0 && start_index < 8);
  assert(stop_index >= 0 && stop_index < 8);
  assert(start_doy < stop_doy);
#endif
  // get average
  int num_ap = 0, curi=start_index;
  double sum_ap = 0e0;
  cur = start_doy;
  while (!(cur == stop_doy && curi == stop_index)) {
#ifdef DEBUG
    printf("+(%ld,%d)=%d", cur - today, curi, cur->ApIndexes[curi]);
#endif
    sum_ap += cur->ApIndexes[curi++];
    if (curi > 7) {
      curi = 0;
      ++cur;
    }
    ++num_ap;
  }
  in.aparr.a[5] = sum_ap / num_ap;
#ifdef DEBUG
  printf("\nAdded %d number of Ap indexes\n", num_ap);
#endif

  // (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
  // TO CURRENT TIME
  // 57 hours = 48 + 9
  // 36 hours = where we last stoped (above)
  stop_index = start_index;
  stop_doy = start_doy;
  start_index = ihours(in.sec) / 3 - 3 + 1;
  start_doy = flux_data_ + 3 - 2; // before previous day
  if (start_index < 0) {
    --start_doy;
    start_index += 7;
  }
#ifdef DEBUG
  printf("\tWill start at %ld days before and index %d\n", start_doy - today, start_index);
  printf("\tWill stop at %ld days before and index %d\n", stop_doy - today, stop_index);
  assert(start_index >= 0 && start_index < 8);
  assert(stop_index >= 0 && stop_index < 8);
  assert(start_doy < stop_doy);
#endif
  // get average
  num_ap = 0;
  sum_ap = 0e0;
  while (!(start_doy == stop_doy && start_index == stop_index)) {
#ifdef DEBUG
    printf("+(%ld,%d)=%d", start_doy - today, start_index, start_doy->ApIndexes[start_index]);
#endif
    sum_ap += start_doy->ApIndexes[start_index++];
    if (start_index > 7) {
      ++start_doy;
      start_index = 0;
    }
    ++num_ap;
  }
  in.aparr.a[6] = (double)sum_ap / (double)num_ap;
#ifdef DEBUG
  printf("\nAdded %d number of Ap indexes\n", num_ap);
#endif

  chi = ihours(in.sec) / 3;
  return 0;
}

int dso::nrlmsise00::detail::Nrlmsise00DataFeed::update(
    dso::nrlmsise00::detail::InParamsCore &in) noexcept {

  dso::modified_julian_day inmjd (dso::ydoy2mjd(in.year, in.doy));
  if (inmjd == mjd_) {
    // date not changed, here comes the data. let's see if the hour interval
    // changed
    if (ihours(in.sec) / 3 == chi) {
      // perfect! nothing to do!
      return 0;
    } else {
      // shit, hour index changed
      printf("-------> shit, computing new indexes; last=%d, new=%d\n", chi, ihours(in.sec) / 3);
#ifdef DEBUG
      assert(chi != 7);
#endif
      this->init(in);
    }
  } else {
    // date changed ! TODO
    return 1;
  }

  return 1;
}
