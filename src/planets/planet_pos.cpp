#include "planetpos.hpp"

dso::iStatus dso::planet_pos(dso::Planet p, const dso::TwoPartDate &mjd_tt,
                        Eigen::Matrix<double, 3, 1> &pos) noexcept {
  /* TT to ET (CSPICE internal system) */
  const double et = dso::cspice::mjdtt2et(mjd_tt);
  
  /* planet's NAIF id */
  int id=0;
  if (dso::cspice::planet_to_naif_id(p, id)) {
    fprintf(
        stderr,
        "[ERROR] Failed to match given plnet to a NAIF id (traceback: %s)\n",
        __func__);
    return dso::iStatus(1);
  }

  double kmpos[3]; /* planet pos in [km] */
  dso::cspice::j2planet_pos_from(et, id, 399 /* i.e. Earth */, kmpos);

  /* assign to output matrix */
  pos(0) = kmpos[0] * 1e3;
  pos(1) = kmpos[1] * 1e3;
  pos(2) = kmpos[2] * 1e3;

  /* all done */
  return dso::iStatus(0);
}
