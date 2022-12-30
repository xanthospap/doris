#include "tides.hpp"
#include <cassert>
#include <cstdio>

struct {
  char d[8];
  int error;
  int iar[6];
  bool check_iarray(const int *iarr) const noexcept {
    for (int i = 0; i < 6; i++)
      if (iarr[i] != iar[i])
        return i + 1;
    return 0;
  }
} TestData[] = {
    {"125.755", 0, {1, -3, 0, 2, 0, 0}}, {"135,655", 0, {1, -2, 0, 1, 0, 0}},
    {"195,465", 0, {1, 4, 0, -1, 1, 0}}, {"167,555", 0, {1, 1, 2, 0, 0, 0}},
    {"073,555", 0, {0, 2, -2, 0, 0, 0}}, {"095,355", 0, {0, 4, 0, -2, 0, 0}},
    {"95,355", 1, {0, 4, 0, -2, 0, 0}},  {"X95,355", 1, {0, 4, 0, -2, 0, 0}},
    {"095/355", 1, {0, 4, 0, -2, 0, 0}}, {"095X355", 1, {0, 4, 0, -2, 0, 0}},
    {"095.35X", 1, {0, 4, 0, -2, 0, 0}}, {"095.35a", 1, {0, 4, 0, -2, 0, 0}},
    {"095.35", 1, {0, 4, 0, -2, 0, 0}},
};
constexpr const int N = sizeof(TestData) / sizeof(TestData[0]);

int main() {
  printf("Unit test for Resolving Doodson Numbers\n");

  int error;
  int array[6];
  for (int i = 0; i < N; i++) {
    error = dso::doodson2intarray(TestData[i].d, array);
    if ((bool)error != (bool)TestData[i].error) {
      fprintf(stderr,
              "ERROR Converting Doodson string [%s]; expected %s got %s\n",
              TestData[i].d, (TestData[i].error) ? "error" : "ok",
              (error) ? "error" : "ok");
      assert(false);
    }
    if (!error) {
      if (TestData[i].check_iarray(array)) {
        fprintf(stderr,
                "ERROR Erronuous conversion of Doodson string for [%s]!\n",
                TestData[i].d);
        fprintf(stderr, "      Expected: ");
        for (int j=0; j<6; j++) fprintf(stderr, "%+d", TestData[i].iar[j]);
        fprintf(stderr, "      got: ");
        for (int j=0; j<6; j++) fprintf(stderr, "%+d", array[j]);
        fprintf(stderr, "\n");
        assert(false);
      }
    }
  }

  printf("All tests passed!\n");
  return 0;
}
