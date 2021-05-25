#include <iostream>
#include <cstdint>
#include <limits>
#include <random>
#include <cassert>
#include <chrono>

constexpr double MISSING = std::numeric_limits<double>::min();
constexpr int MAX_UINT8_T = std::numeric_limits<uint8_t>::max();
constexpr int MIN_UINT8_T = std::numeric_limits<uint8_t>::min();
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> distrib(MIN_UINT8_T, MAX_UINT8_T);
std::uniform_real_distribution<> dis(MISSING, std::numeric_limits<double>::max());
static int counter = 0;

struct Measurement1 {
  double val_;
  int_fast8_t missing_, 
    flagm1_, 
    flagm2_;
};
Measurement1 rand_m1() {
  ++counter;
  Measurement1 m;
  m.val_ = dis(gen);
  m.missing_ = 0;
  if (!(counter/33)) m.missing_ = 1;
  m.flagm1_ = distrib(gen);
  m.flagm2_ = distrib(gen);
  return m;
}
double foo1(const Measurement1& m) noexcept {
  if (m.missing_) {
    return -999.99e0;
  } else {
    if (m.flagm1_ > 16) {
      return m.val_ + 1e3;
    } else if (m.flagm2_ < 16) {
      return m.val_ - 1e-3;
    } else {
      return m.val_ / 1e3;
    }
  }
}

struct Measurement2 {
  double val_;
  int_fast8_t flagm1_, 
    flagm2_;
};
double foo1(const Measurement2& m) noexcept {
  if (m.val_ == MISSING) {
    return -999.99e0;
  } else {
    if (m.flagm1_ > 16) {
      return m.val_ + 1e3;
    } else if (m.flagm2_ < 16) {
      return m.val_ - 1e-3;
    } else {
      return m.val_ / 1e3;
    }
  }
}
Measurement2 rand_m2() {
  ++counter;
  Measurement2 m;
  m.val_ = dis(gen);
  if (!(counter/33)) m.val_ = MISSING;
  m.flagm1_ = distrib(gen);
  m.flagm2_ = distrib(gen);
  return m;
}

struct Flag {
  uint_fast32_t flag_{0};
  // flag1, flag2, missing
  uint_fast8_t get_flag1() const noexcept {
    return flag_ & 0xFF;
  }
  uint_fast8_t get_flag2() const noexcept {
    return (flag_ >> 8) & 0xFF;
  }
  uint_fast8_t get_missing() const noexcept {
    return (flag_ >> 16) & 0xFF;
  }
  void set_flag1(uint8_t val) noexcept {
    uint_fast8_t a(0x0), b(get_missing()),c(get_flag2());
    flag_ = (a << 24) | (b << 16) | (c << 8) | val;
  }
  void set_flag2(uint8_t val) noexcept {
    uint_fast8_t a(0x0), b(get_missing()), c(val), d(get_flag1());
    flag_ = (a << 24) | (b << 16) | (c << 8) | d;
  }
  void set_missing(uint8_t val) noexcept {
    uint_fast8_t a(0x0), b(val), c(get_flag2()), d(get_flag1());
    flag_ = (a << 24) | (b << 16) | (c << 8) | d;
  }
};
struct Measurement3 {
  double val_;
  Flag flag_;
};
double foo1(const Measurement3& m) noexcept {
  if (m.flag_.get_missing()) {
    return -999.99e0;
  } else {
    if (m.flag_.get_flag1() > 16) {
      return m.val_ + 1e3;
    } else if (m.flag_.get_flag2() < 16) {
      return m.val_ - 1e-3;
    } else {
      return m.val_ / 1e3;
    }
  }
}
Measurement3 rand_m3() {
  ++counter;
  Measurement3 m;
  m.val_ = dis(gen);
  if (!(counter/33)) m.flag_.set_missing(1);
  m.flag_.set_flag1(distrib(gen));
  m.flag_.set_flag2(distrib(gen));
  return m;
}

using namespace std::chrono;
int main() {
  printf("Here are the struct sizes:\n");
  printf("Size of Measurement1:  %u\n", sizeof(Measurement1));
  printf("Size of Measurement2:  %u\n", sizeof(Measurement2));
  printf("Size of Measurement3:  %u\n", sizeof(Measurement3));

  Flag f;
  printf("Trying out the Flag struct\n");
  printf("set flag1 to 9\n");
  f.set_flag1(9);
  printf("set flag2 to 23\n");
  f.set_flag2(23);
  printf("set missing to 1\n");
  f.set_missing(1);
  printf("flag1   -> %d\n", f.get_flag1());
  printf("flag2   -> %d\n", f.get_flag2());
  printf("missing -> %d\n", f.get_missing());
  printf("Let's do that one million times to be sure ...");
  for (int i=0; i<1000000; ++i) {
    uint_fast8_t f1=distrib(gen), f2=distrib(gen), msng=distrib(gen);
    f.set_flag1(f1);
    f.set_flag2(f2);
    f.set_missing(msng);
    assert(f.get_flag1() == f1);
    assert(f.get_flag2() == f2);
    assert(f.get_missing() == msng);
  }
  printf(" cool, everything looks ok!\n");

  for (int j=0; j<10; ++j) {
    auto start = high_resolution_clock::now();
    for (int i=0; i<1000000; ++i) {
      auto m = rand_m1();
      auto d = foo1(m);
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken by version 1: " << duration.count() << " microseconds\n";

    counter = 0;
    start = high_resolution_clock::now();
    for (int i=0; i<1000000; ++i) {
      auto m = rand_m2();
      auto d = foo1(m);
    }
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken by version 2: " << duration.count() << " microseconds\n";

    counter = 0;
    start = high_resolution_clock::now();
    for (int i=0; i<1000000; ++i) {
      auto m = rand_m3();
      auto d = foo1(m);
    }
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken by version 3: " << duration.count() << " microseconds\n";
  }

  return 0;
}
