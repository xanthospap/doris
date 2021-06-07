
void cartesian2keplerian_tapley(const double *x, const double *xdot) noexcept {

  const double rsq = dot_product(x, x);
  const double r = std::sqrt(rsq);
  double hv[3];
  cross_product(x, xdot, hv);
  const double hsq = dot_product(hv, hv);
  const double h = std::sqrt(hsq);

  double i = std::acos(hv[2] / h);
  double Omega = std::atan2(hv[0], -hv[1]);

  double ksi = dot_product(xdot, xdot)/2e0 - ngpt::GM / r;
  double a = -ngpt::GM / ksi / 2e0;
  double e = std::sqrt(1e0 + (2e0*ksi*hsq) / ngpt::GM / ngpt::GM);
  double p = hsq / ngpt::GM;

  double cosf = (p-r) / r / e;
  double sinf = (p/h/e) * dot_product(x, xdot) / r;
  double cosof = (x[0]/r)*std::cos(Omega) + (x[1]/r)*std::sin(Omega);
  double sinof = x[2] / r / std::sin(i);
  double fangle = std::atan(sinf / cosf);
  double ofangle = std::atan(sinof / cosof);
  double omega = ofangle - fangle;

  double cosE = r * cosf / a + e;
  double b = a * std::sqrt(1e0 - e*e);
  double sinE = r * sinf / b;
  double E = std::atan(sinE / cosE);
  double M = E - e * sinE;
  
  printf("a    = %15.6f\n", a);
  printf("e    = %15.13f\n", e);
  printf("ω    = %15.12f\n", ngpt::rad2deg(omega));
  printf("Ω    = %15.12f\n", ngpt::rad2deg(Omega));
  printf("i    = %15.12f\n", ngpt::rad2deg(i));
  // printf("T0   = %15.12f\n", T0);
  printf("M    = %15.12f\n", ngpt::rad2deg(M));
}

void cartesian2keplerian_beutler(const double *x, const double *xdot,
                                 double t = 0e0) noexcept {
  /*
   * Output:
   * -----------------------------------
   * semi-latus rectum, p
   * eccentricity, e
   * argument of periapsis, ω (rad)
   * longitude of ascending node, Ω (rad)
   * inclination, i (rad)
   * time of perihelion, T0
   */
  double hv[3];
  cross_product(x, xdot, hv);
  hv[1] = -hv[1]; // bernese-like cross product
  double hsq = dot_product(hv, hv);
  double h = std::sqrt(hsq);

  // ascending node, Ω
  double Omega = std::atan2(hv[0], hv[1]);

  // inclination, i
  double i = std::acos(hv[2] / h);

  // constant Edot
  double rsq = dot_product(x, x);
  double r = std::sqrt(rsq);
  double Edot = dot_product(xdot, xdot) / 2e0 - ngpt::GM / r;

  // argument p
  double p = hsq / ngpt::GM;

  // argument e
  double e = std::sqrt(1e0 + 2e0 * hsq * Edot / ngpt::GM / ngpt::GM);

  // compute the argument of latitude at time t
  // this is: R_x(i) * R_z(Ω) * r (expanded to ...)
  double ck = std::cos(Omega);
  double sk = std::sin(Omega);
  double ci = std::cos(i);
  double si = std::sin(i);
  double xx0 = ck * x[0] + sk * x[1];
  double xx1 = -ci * sk * x[0] + ci * ck * x[1] + si * x[2];
  // double xx2 = si*sk*x[0]-si*ck*x[1]+ci*x[2];
  double u = std::atan2(xx1, xx0);
  double ecv = p / r - 1e0;
  // double esv = h / ngpt::GM;
  double esv = std::sqrt(p / ngpt::GM) / r * dot_product(x, xdot);
  double v1 = std::atan2(esv, ecv);
  double omega = u - v1;
  // for elliptic orbits (eq. 5.33)
  double a = p / (1e0 - e * e);
  double sinv = std::sin(v1);
  double cosv = std::cos(v1);
  /*
  double sinE = (std::sqrt(1e0-e*e) / (1e0+e)) * (sinv/cosv); // 4.51
  double cosE = (e + cosv) / (1e0+e*cosv);                    // 4.51
  double E = std::atan(sinE/cosE);
  */
  double E = std::atan(std::sqrt(1e0 + e * e) * sinv / (e + cosv));
  double T0 =
      t - (E - e * std::sin(E)) / std::sqrt(ngpt::GM / a / a / a); // Table 4.2

  printf("a    = %15.6f\n", a);
  printf("e    = %15.13f\n", e);
  printf("ω    = %15.12f\n", ngpt::rad2deg(omega));
  printf("Ω    = %15.12f\n", ngpt::rad2deg(Omega));
  printf("i    = %15.12f\n", ngpt::rad2deg(i));
  printf("T0   = %15.12f\n", T0);
  // printf("M    = %15.12f\n", ngpt::rad2deg(M));
}

void cartesian2keplerian_bernese(const double *x, const double *v,
                                 double t = 0e0) noexcept {
  /*
   * Output:
   * -----------------------------------
   * semi-major axis, a (m)
   * eccentricity, e
   * argument of periapsis, ω (rad) -- in bernese PER --
   * longitude of ascending node, Ω (rad) -- in bernese KN --
   * inclination, i (rad)
   * time of perigee/helion passing, T0
   */
  double h[3] = {
      x[1] * v[2] - x[2] * v[1],
      -x[2] * v[0] + x[0] * v[2],
      x[0] * v[1] - x[1] * v[0],
  };

  double kn = std::atan2(h[0], h[1]);

  double i = std::atan2(std::sqrt(h[0] * h[0] + h[1] * h[1]), h[2]);

  double p = (h[0] * h[0] + h[1] * h[1] + h[2] * h[2]) / ngpt::GM;
  double r = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  double ecv = p / r - 1e0;
  double esv =
      std::sqrt(p / ngpt::GM) / r * (x[0] * v[0] + x[1] * v[1] + x[2] * v[2]);
  double v1 = std::atan2(esv, ecv);
  double e = std::sqrt(ecv * ecv + esv * esv);

  double ck = std::cos(kn);
  double sk = std::sin(kn);
  double ci = std::cos(i);
  double si = std::sin(i);
  double xx[] = {ck * x[0] + sk * x[1],
                 -ci * sk * x[0] + ci * ck * x[1] + si * x[2],
                 si * sk * x[0] - si * ck * x[1] + ci * x[2]};

  double u = std::atan2(xx[1], xx[0]);
  double per = u - v1;
  double ex = 2e0 * std::atan(std::sqrt((1e0 - e) / (1e0 + e)) *
                              (std::sin(v1 / 2) / std::cos(v1 / 2e0)));
  double a = p / (1e0 - e * e);
  double a3 = a * a * a;
  double t0 = t - (ex - e * std::sin(ex)) / std::sqrt(ngpt::GM / a3);

  printf("a    = %15.6f\n", a);
  printf("e    = %15.13f\n", e);
  printf("ω    = %15.12f\n", ngpt::rad2deg(per));
  printf("Ω    = %15.12f\n", ngpt::rad2deg(kn));
  printf("i    = %15.12f\n", ngpt::rad2deg(i));
  printf("T0   = %15.12f\n", t0);
  // printf("M    = %15.12f\n", ngpt::rad2deg(M));
  return;
}

/// https://downloads.rene-schwarz.com/download/M002-Cartesian_State_Vectors_to_Keplerian_Orbit_Elements.pdf
void cartesian2keplerian_schwarz(const double *x, const double *v) noexcept {
  /*
   * Output:
   * -----------------------------------
   * semi-major axis, a (m)
   * eccentricity, e
   * argument of periapsis, ω (rad)
   * longitude of ascending node, Ω (rad)
   * inclination, i (rad)
   * mean anomaly, M (rad)
   */

  // Calculate orbital momentum vector
  double hv[3];
  cross_product(x, v, hv);
  const double hsq = dot_product(hv, hv);
  const double h = std::sqrt(hsq);
  const double rsq = dot_product(x, x);
  const double r = std::sqrt(rsq);

  // eccentricity vector
  double ev[3];
  cross_product(v, hv, ev);
  ev[0] /= ngpt::GM;
  ev[1] /= ngpt::GM;
  ev[2] /= ngpt::GM;
  ev[0] -= (x[0] / r);
  ev[1] -= (x[1] / r);
  ev[2] -= (x[2] / r);
  const double esq = dot_product(ev, ev);

  // eccentricity, e
  double e = std::sqrt(esq);

  // Determine the vector n pointing towards the ascending node
  double nv[] = {-hv[1], hv[0], 0e0};
  double n = std::sqrt(dot_product(nv, nv));

  //  True anomaly
  double ta = std::acos(dot_product(ev, x) / e / r);
  if (dot_product(x, v) < 0e0)
    ta = ngpt::D2PI - ta;

  // Orbital inclination, i
  double i = std::acos(hv[2] / h);

  // Eccentric anomaly, E
  double E =
      2e0 * std::atan2(std::tan(ta / 2e0), std::sqrt((1e0 + e) / (1e0 - e)));
  // if (E<0e0) E += ngpt::D2PI;

  // longtitude of the acending node, Ω
  double Omega = std::acos(nv[0] / n);
  // if (nv[1]<0e0) Omega = ngpt::D2PI - Omega;

  // argument of periapsis, ω
  double omega = std::acos(dot_product(nv, ev) / n / e);
  // if (ev[2]<0e0) omega = ngpt::D2PI - omega;

  // Mean anomaly, M
  double M = E - e * std::sin(E);

  // Semi-major axis, a
  double a = 1e0 / ((2e0 / r) - (dot_product(v, v) / ngpt::GM));

  printf("a    = %15.6f\n", a);
  printf("e    = %15.13f\n", e);
  printf("ω    = %15.12f\n", ngpt::rad2deg(omega));
  printf("Ω    = %15.12f\n", ngpt::rad2deg(Omega));
  printf("i    = %15.12f\n", ngpt::rad2deg(i));
  printf("M    = %15.12f\n", ngpt::rad2deg(M));
}

int main() {
  double x[] = {5492000.34e0, 3984001.40e0, 2955.81e0};
  double v[] = {-3931.046491e0, 5498.676921, 3665.980697};

  printf("/*--Schwarz "
         "Impl.----------------------------------------------------------------"
         "-*/\n");
  cartesian2keplerian_schwarz(x, v);
  printf("/*--Beutler.---------------------------------------------------------"
         "-------------*/\n");
  cartesian2keplerian_beutler(x, v);
  printf("/*--Bernese.---------------------------------------------------------"
         "-------------*/\n");
  cartesian2keplerian_bernese(x, v);
  printf("/*--Tapley .---------------------------------------------------------"
         "-------------*/\n");
  cartesian2keplerian_tapley(x, v);

  double x2[] = {-6045000.0e0, -3490000.0e0, +2500000.0e0};
  double v2[] = {-3457.0e0, +6618.0e0, 2533.0e0};

  printf("/*--Schwarz "
         "Impl.----------------------------------------------------------------"
         "-*/\n");
  cartesian2keplerian_schwarz(x2, v2);
  printf("/*--Beutler.---------------------------------------------------------"
         "-------------*/\n");
  cartesian2keplerian_beutler(x2, v2);
  printf("/*--Bernese.---------------------------------------------------------"
         "-------------*/\n");
  cartesian2keplerian_bernese(x2, v2);
  printf("/*--Tapley .---------------------------------------------------------"
         "-------------*/\n");
  cartesian2keplerian_tapley(x2, v2);
  return 0;
}
