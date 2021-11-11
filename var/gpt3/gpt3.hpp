#ifndef __GPT3_TROPO_HPP__
#define __GPT3_TROPO_HPP__

namespace dso {
    namespace gpt3 {
constexpr int GPT3_5_GRID_LINES = 2592;

struct gpt3_5_grid {
  double p_grid[GPT3_5_GRID_LINES][5];
  double T_grid[GPT3_5_GRID_LINES][5];
  double Q_grid[GPT3_5_GRID_LINES][5];
  double dT_grid[GPT3_5_GRID_LINES][5];
  double u_grid[GPT3_5_GRID_LINES];
  double Hs_grid[GPT3_5_GRID_LINES];
  double ah_grid[GPT3_5_GRID_LINES][5];
  double aw_grid[GPT3_5_GRID_LINES][5];
  double la_grid[GPT3_5_GRID_LINES][5];
  double Tm_grid[GPT3_5_GRID_LINES][5];
  double Gn_h_grid[GPT3_5_GRID_LINES][5];
  double Ge_h_grid[GPT3_5_GRID_LINES][5];
  double Gn_w_grid[GPT3_5_GRID_LINES][5];
  double Ge_w_grid[GPT3_5_GRID_LINES][5];

  void ugrid_slice(const int *indexes, double *values) noexcept {
      values[0] = u_grid[indexes[0]];
      values[1] = u_grid[indexes[1]];
      values[2] = u_grid[indexes[2]];
      values[3] = u_grid[indexes[3]];
      return;
  }
  void hsgrid_slice(const int *indexes, double *values) noexcept {
      values[0] = Hs_grid[indexes[0]];
      values[1] = Hs_grid[indexes[1]];
      values[2] = Hs_grid[indexes[2]];
      values[3] = Hs_grid[indexes[3]];
      return;
  }


}; // gpt3_5_grid

int parse_gpt3_5_grid(const char *gridfn, gpt3_5_grid *grid) noexcept;
    }// gpt3
}// dso

#endif