#include <cassert>
#include <cstdio>

struct BaseInterpolant {
  const double *x, *y;
  int size;
  BaseInterpolant(const double *xin, const double *yin, int sz) noexcept
      : x(xin), y(yin), size(sz){};
};

struct NearestNeighborInterpolator {
  const BaseInterpolant base;
  NearestNeighborInterpolator(const double *xin, const double *yin,
                              int sz) noexcept
      : base(xin, yin, sz){};
  double value_at(double x, int& status) const noexcept {
    status = 0;
    if (x <= base.x[0]) {
      if (x < base.x[0]) {
        status = 1;
      }
      return base.y[0];
    }
    for (int i = 1; i < base.size; ++i) {
      if (x >= base.x[i - 1] && x < base.x[i]) {
        double dist_xt1 = x - base.x[i - 1];
        double dist_xt2 = base.x[i] - x;
        return (dist_xt1 < dist_xt2) ? base.y[i - 1] : base.y[i];
      }
    }
    status = 2;
    return base.y[base.size - 1];
  }
};

struct LinearInterpolator {
  const BaseInterpolant base;
  LinearInterpolator(const double *xin, const double *yin, int sz) noexcept
      : base(xin, yin, sz){};
  double value_at(double x, int& status) const noexcept {
    status = 0;
    if (x <= base.x[0]) {
      if (x < base.x[0]) {
        status = 1;
      }
      return base.y[0];
    }
    for (int i = 1; i < base.size; ++i) {
      if (x >= base.x[i - 1] && x < base.x[i]) {
        return base.y[i - 1] +
               ((base.y[i - 1] - base.y[i]) / (base.x[i - 1] - base.x[i])) *
                   (x - base.x[i - 1]);
      }
    }
    status = 2;
    return base.y[base.size - 1];
  }
};

struct QuadraticInterpolator {
  const BaseInterpolant base;
  QuadraticInterpolator(const double *xin, const double *yin, int sz) noexcept
      : base(xin, yin, sz){};
  double value_at(double x, int& status) const noexcept {
    status = 0;
    if (x <= base.x[0]) {
      if (x < base.x[0]) {
        status = 1;
      }
      return base.y[0];
    }
    for (int i = 1; i < base.size; ++i) {
      if (x >= base.x[i - 1] && x < base.x[i]) {
        int idx = -1;
        if (i > 1 && i < base.size - 1) {
          double dist_xt1 = x - base.x[i - 2];
          double dist_xt2 = base.x[i + 1] - x;
          idx = (dist_xt1 < dist_xt2) ? i - 2 : i + 1;
        } else if (i <= 1 && i < base.size - 1) {
          idx = i + 1;
        } else if (i > 1 && i >= base.size - 1) {
          idx = i - 2;
        }
        assert(idx >= 0 && idx <= base.size - 1);
        // printf("\n\tInterpolation points: [%1d, %1d, %1d]\n", i-1, i, idx);
        double y0 = base.y[i - 1];
        double y1 = base.y[i];
        double x0 = base.x[i - 1];
        double x1 = base.x[i];
        double y2 = base.y[idx];
        double x2 = base.x[idx];
        double L0 = ((x - x1) * (x - x2)) / ((x0 - x1) * (x0 - x2));
        double L1 = ((x - x0) * (x - x2)) / ((x1 - x0) * (x1 - x2));
        double L2 = ((x - x0) * (x - x1)) / ((x2 - x0) * (x2 - x1));
        return y0 * L0 + y1 * L1 + y2 *L2;
      }
    }
    status = 2;
    return base.y[base.size - 1];
  }
};

int main() {
  constexpr int size = 90 / 5 + 1;
  double xpts[size] = {0},
         ypts[size] = {0.00e0,  2.05e0,  7.24e0,  9.21e0,  6.71e0,
                       8.14e0,  11.87e0, 12.48e0, 12.28e0, 13.67e0,
                       13.91e0, 13.01e0, 13.01e0, 11.87e0, 9.70e0,
                       7.94e0,  4.99e0,  0.41e0,  -3.93};
  for (int i = 1; i < size; i++) xpts[i] = xpts[i - 1] + 5e0;

  for (int i = 0; i < size; i++) printf("%15.5f %15.5f\n", xpts[i], ypts[i]);

  constexpr int test_size = 160;
  constexpr double step = 100 / (double)160;
  constexpr double start_x = -0.5e0;
  // double new_x[test_size], new_y[test_size];
  // int new_status[test_size];

  NearestNeighborInterpolator nintrp(xpts, ypts, size);
  LinearInterpolator lintrp(xpts, ypts, size);
  QuadraticInterpolator qintrp(xpts, ypts, size);
  for (int i = 0; i < test_size; i++) {
    int nstat, lstat, qstat;
    double x = start_x + step * i;
    double yn = nintrp.value_at(x, nstat);
    double yl = lintrp.value_at(x, lstat);
    double yq = qintrp.value_at(x, qstat);
    printf("%15.5f %15.5f %1d %15.5f %1d %15.5f %1d\n", x, yn, nstat, yl, lstat, yq, qstat);
  }

  return 0;
}
