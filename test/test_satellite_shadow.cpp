#include "datetime/datetime_write.hpp"
#include "planetpos.hpp"
#include "satellite.hpp"
#include "sp3/sp3.hpp"
#include "sp3/sv_interpolate.hpp"
#include <cassert>
#include <cstdio>
#include <algorithm>

bool transform2inertial = false;

constexpr const dso::datetime<dso::nanoseconds> tmax =
    dso::datetime<dso::nanoseconds>::max();

char buf[64];
char *pd(const dso::datetime<dso::nanoseconds> &t) {
  if (t == tmax)
    std::strcpy(buf, "---");
  else
    std::strcpy(buf, dso::strftime_ymd_hmfs(t).c_str());
  return buf;
}

struct ShadowPass {
  ShadowPass(dso::datetime<dso::nanoseconds> t = tmax)
      : pstart{t}, pstop{tmax}, ustart{tmax}, ustop{tmax} {};
  dso::datetime<dso::nanoseconds> pstart{tmax};
  dso::datetime<dso::nanoseconds> pstop{tmax}; // penumbra
  dso::datetime<dso::nanoseconds> ustart{tmax};
  dso::datetime<dso::nanoseconds> ustop{tmax}; // umbra
  bool in_umbra() const noexcept { return (ustart != tmax && ustop == tmax); }
  bool in_penumbra() const noexcept {
    return (pstart != tmax && pstop == tmax);
  }
  bool arc_closed() const noexcept { return pstop != tmax; }
  bool arc_is_empty() const noexcept {
    return (pstart == tmax && pstop == tmax && ustart == tmax && ustop == tmax);
  }
  void print() const {
    dso::nanoseconds dns = dso::delta_sec(pstop, pstart);
    printf("%s %25s %25s %s %25s %25s %.3f\n", "Penumbra", pd(pstart),
           pd(pstop), "Umbra", pd(ustart), pd(ustop),
           dns.to_fractional_seconds());
  }
};

void shadow_coeff(const dso::Sp3DataBlock *block, const double *vsun) noexcept {
  dso::Vector3 rsat{
      {block->state[0] * 1e3, block->state[1] * 1e3, block->state[2] * 1e3}};
  dso::Vector3 rsun{{vsun[0] * 1e3, vsun[1] * 1e3, vsun[2] * 1e3}};
  double c1 = utest::montebruck_shadow(rsat, rsun);
  double c2 = utest::conical_shadow(rsat, rsun);
  double c3 = utest::vallado_shadow(rsat, rsun);
  printf("%.5f f1=%.3f f2=%.3f f3=%.3f\n", block->t.as_mjd(), c1, c2, c3);
}

void shadow_coeff(const dso::datetime<dso::nanoseconds> &t,
                  const double *sv_xyz, const double *vsun,
                  std::vector<ShadowPass> &montebruck,
                  std::vector<ShadowPass> &vallado,
                  std::vector<ShadowPass> &conic,
                  std::vector<ShadowPass> &bernese1) noexcept {
  static std::size_t count = 0;
  if (count == 0)
    assert(vallado[0].pstart == tmax);
  dso::Vector3 rsat{{sv_xyz[0] * 1e3, sv_xyz[1] * 1e3, sv_xyz[2] * 1e3}};
  dso::Vector3 rsun{{vsun[0] * 1e3, vsun[1] * 1e3, vsun[2] * 1e3}};

  double c1 = utest::montebruck_shadow(rsat, rsun);
  double c2 = utest::conical_shadow(rsat, rsun);
  double c3 = utest::vallado_shadow(rsat, rsun);
  double c4 = utest::bernese_shadow1(rsat, rsun);
  printf("%6lu %.5f f1=%.3f f2=%.3f f3=%.3f f4=%.3f\n", count, t.as_mjd(), c1, c2, c3, c4);

  auto idx = montebruck.size();
  // if c1 == 1, we are not in the shadow; no discrimination of umbra/penumbra
  if (c1 == 0) { // in penumbra ...
    // new shadow pass ...
    if (montebruck[idx - 1].arc_closed())
      montebruck.emplace_back(t);
  } else { // in light ...
    if (!(montebruck[idx - 1].arc_closed()))
      montebruck[idx - 1].pstop = t;
  }
  
  idx = bernese1.size();
  // if c4 == 1, we are not in the shadow; no discrimination of umbra/penumbra
  if (c4 == 0) { // in penumbra ...
    // new shadow pass ...
    if (bernese1[idx - 1].arc_closed())
      bernese1.emplace_back(t);
  } else { // in light ...
    if (!(bernese1[idx - 1].arc_closed()))
      bernese1[idx - 1].pstop = t;
  }

  idx = vallado.size() - 1;
  if (c3 < 1e0) { // in penumbra or umbra ...
    // printf(">>vallado in shadow ... %.5f\n", c3);
    if (c3 == 0e0) { // in umbra
      // printf("\tvallado in umbra ...\n");
      // are we starting a new arc ?
      if (vallado[idx].arc_is_empty()) {
        // printf("\tlast arc was empty; filling ...\n");
        vallado[idx].pstart = vallado[idx].ustart = t;
      } else if (vallado[idx].pstop != tmax && vallado[idx].pstop != tmax) {
        // printf("\tlast arc was closed; opening new umbra/penubra arc ...\n");
        vallado.emplace_back(t);
        ++idx;
        vallado[idx].ustart = t;
      } else if (vallado[idx].pstart != tmax && vallado[idx].ustart == tmax) {
        // printf("\tlast arc was penumbra; opening new umbra arc ...\n");
        vallado[idx].ustart = t;
      }
    }      // in umbra
    else { // in penumbra
      // printf("\tvallado in penumbra ...\n");
      if (vallado[idx].arc_is_empty()) {
        // printf("\tlast arc was empty; filling ...\n");
        vallado[idx].pstart = t;
      } else if (vallado[idx].pstop == tmax) {
        // printf("\tlast arc was closed; opening new penubra arc ...\n");
        vallado.emplace_back(t);
        ++idx;
      }
    }
  } // in penumbra or umbra
  else {
    if (vallado[idx].pstart != tmax && vallado[idx].pstop == tmax) {
      // printf("\tlast arc was open; closing umbra/penubra arc ...\n");
      vallado[idx].pstop = t;
      if (vallado[idx].ustop == tmax)
        vallado[idx].ustop = t;
    }
  }

  idx = conic.size() - 1;
  if (c3 < 1e0) { // in penumbra or umbra ...
    // printf(">>conic in shadow ... %.5f\n", c3);
    if (c3 == 0e0) { // in umbra
      // printf("\tconic in umbra ...\n");
      // are we starting a new arc ?
      if (conic[idx].arc_is_empty()) {
        // printf("\tlast arc was empty; filling ...\n");
        conic[idx].pstart = conic[idx].ustart = t;
      } else if (conic[idx].pstop != tmax && conic[idx].pstop != tmax) {
        // printf("\tlast arc was closed; opening new umbra/penubra arc ...\n");
        conic.emplace_back(t);
        ++idx;
        conic[idx].ustart = t;
      } else if (conic[idx].pstart != tmax && conic[idx].ustart == tmax) {
        // printf("\tlast arc was penumbra; opening new umbra arc ...\n");
        conic[idx].ustart = t;
      }
    }      // in umbra
    else { // in penumbra
      // printf("\tconic in penumbra ...\n");
      if (conic[idx].arc_is_empty()) {
        // printf("\tlast arc was empty; filling ...\n");
        conic[idx].pstart = t;
      } else if (conic[idx].pstop == tmax) {
        // printf("\tlast arc was closed; opening new penubra arc ...\n");
        conic.emplace_back(t);
        ++idx;
      }
    }
  } // in penumbra or umbra
  else {
    if (conic[idx].pstart != tmax && conic[idx].pstop == tmax) {
      // printf("\tlast arc was open; closing umbra/penubra arc ...\n");
      conic[idx].pstop = t;
      if (conic[idx].ustop == tmax)
        conic[idx].ustop = t;
    }
  }

  ++count;
  return;
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr, "Usage: %s [SP3] [DE/SPK KERNEL] [LEAPSEC/LSK KERNEL]\n",
            argv[0]);
    return 1;
  }

  dso::Sp3c sp3(argv[1]);
#ifdef DEBUG
  sp3.print_members();
  printf("\n");
#endif
  dso::sp3::SatelliteId sv;
  if (sp3.num_sats() != 1) {
    fprintf(stderr, "More than one satellites found in SP3 file! Listing ...\n");
    for (auto i=sp3.sattellite_vector().cbegin(); i<sp3.sattellite_vector().cend(); i++) {
      printf("%3s ", i->to_string().c_str());
    }
    printf("\nChoose one: ");
    char svbuf[64] = {'\0'};
    fgets(svbuf,64,stdin);
    sv.set_id(svbuf);
    if (std::find(sp3.sattellite_vector().cbegin(), sp3.sattellite_vector().cend(), sv) == sp3.sattellite_vector().cend()) {
      fprintf(stderr, "Erronuous satellite id [%s]! not listed in Sp3 file\n", svbuf);
      return 1;
    }
  } else {
    sv.set_id(sp3.sattellite_vector()[0].id);
  }
  // dso::Sp3DataBlock block;

  // setup the interpolation
  dso::SvInterpolator sv_intrp(sv, sp3);
  auto start_t = sp3.start_epoch();
  auto every_t = dso::datetime_interval<dso::nanoseconds>(
      dso::modified_julian_day(0),
      dso::nanoseconds(10 * dso::nanoseconds::sec_factor<long>()));
  const dso::datetime<dso::nanoseconds> *stop_t = sv_intrp.last_block_date();

  // to compute the planet's position via cspice, we need to load:
  // 1. the planetary ephemeris (SPK) kernel
  // 2. the leap-second (aka LSK) kernel
  dso::cspice::load_if_unloaded_spk(argv[2]);
  dso::cspice::load_if_unloaded_lsk(argv[3]);

  std::vector<ShadowPass> sp_mon(1, ShadowPass());
  std::vector<ShadowPass> sp_val(1, ShadowPass());
  std::vector<ShadowPass> sp_con(1, ShadowPass());
  std::vector<ShadowPass> sp_bn1(1, ShadowPass());
  double vsun[3], xyz[3], dxdydz[3];
  assert(sp_val[0].pstart == tmax);
  while (start_t < *stop_t) {
    // interpolate SV position at time t
    if (!sv_intrp.interpolate_at(start_t, xyz, dxdydz)) {

      // sun position vector, geocentric [km]
      dso::sun_vector_cspice(start_t, vsun);

      //if (transform2inertial) {
      //  auto ut1 = t;
      //  ut1.add_seconds(dso::nanoseconds(69 * 1'000'000'000L));
      //  auto c2t = iers2010::sofa::c2t06a(t.as);
      //}

      // compute shadow functions ...
      shadow_coeff(start_t, xyz, vsun, sp_mon, sp_val, sp_con, sp_bn1);
    }

    start_t += every_t;
  }
  printf("Note, here are the sat coordinates: %.3f %.3f %.3f\n", xyz[0], xyz[1],
         xyz[2]);

  // print Shadow intervals ...
  printf("%5s %25s %25s\n", "Arc", "Penumbra", "Umbra");
  printf("Shadow Interval/Passes for model Montenbruck:\n");
  for (auto it = sp_mon.cbegin(); it != sp_mon.cend(); ++it) {
    // printf("%5d %.5f %.5f\n", (int)std::distance(sp_mon.cbegin(), it),
    // it->pstart.as_mjd(), it->pstop.as_mjd());
    printf("%5d ", (int)std::distance(sp_mon.cbegin(), it));
    it->print();
  }
  printf("Shadow Interval/Passes for model Vallado:\n");
  for (auto it = sp_val.cbegin(); it != sp_val.cend(); ++it) {
    // printf("%5d %.5f %.5f %.5f %.5f\n", (int)std::distance(sp_val.cbegin(),
    // it), it->pstart.as_mjd(), it->pstop.as_mjd(), it->ustart.as_mjd(),
    // it->ustop.as_mjd());
    printf("%5d ", (int)std::distance(sp_val.cbegin(), it));
    it->print();
  }
  printf("Shadow Interval/Passes for model Conic:\n");
  for (auto it = sp_con.cbegin(); it != sp_con.cend(); ++it) {
    // printf("%5d %.5f %.5f %.5f %.5f\n", (int)std::distance(sp_con.cbegin(),
    // it), it->pstart.as_mjd(), it->pstop.as_mjd(), it->ustart.as_mjd(),
    // it->ustop.as_mjd());
    printf("%5d ", (int)std::distance(sp_con.cbegin(), it));
    it->print();
  }
  printf("Shadow Interval/Passes for model Bernese-1:\n");
  for (auto it = sp_bn1.cbegin(); it != sp_bn1.cend(); ++it) {
    // printf("%5d %.5f %.5f %.5f %.5f\n", (int)std::distance(sp_con.cbegin(),
    // it), it->pstart.as_mjd(), it->pstop.as_mjd(), it->ustart.as_mjd(),
    // it->ustop.as_mjd());
    printf("%5d ", (int)std::distance(sp_bn1.cbegin(), it));
    it->print();
  }

  printf("number of shadow passes for odel montenbruck: %d\n", (int)sp_mon.size());
  printf("number of shadow passes for odel vallado    : %d\n", (int)sp_val.size());
  printf("number of shadow passes for odel conic      : %d\n", (int)sp_con.size());
  printf("number of shadow passes for odel bernese-1  : %d\n", (int)sp_bn1.size());

  // let's try reading the records; note that -1 denotes EOF
  // int j;
  // std::size_t rec_count = 0;
  // do {
  //  j = sp3.get_next_data_block(sv, block);
  //  if (j > 0) {
  //    fprintf(stderr, "Something went wrong ....status = %3d\n", j);
  //    return 1;
  //  } else if (j == -1) {
  //    printf("EOF encountered; Sp3 file read through!\n");
  //  }

  //  // sun position vector, geocentric [km]
  //  dso::sun_vector_cspice(block.t, vsun);
  //
  //  // compute shadow functions ...
  //  shadow_coeff(&block, vsun);

  //  ++rec_count;
  //} while (!j);

  // printf("Num of records read: %6lu\n", rec_count);

  return 0;
}
