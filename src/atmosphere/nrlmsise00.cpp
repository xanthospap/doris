#include <cmath>
#include <cstdio>

/* PARMB */
static double gsurf;
static double re;

/* GTS3C */
static double dd;

/* DMIX */
static double dm04, dm16, dm28, dm32, dm40, dm01, dm14;

/* MESO7 */
static double meso_tn1[5];
static double meso_tn2[4];
static double meso_tn3[5];
static double meso_tgn1[2];
static double meso_tgn2[2];
static double meso_tgn3[2];

/* POWER7 */
extern double pt[150];
extern double pd[9][150];
extern double ps[150];
extern double pdl[2][25];
extern double ptl[4][100];
extern double pma[10][100];
extern double sam[100];

/* LOWER7 */
extern double ptm[10];
extern double pdm[8][10];
extern double pavgm[10];

/* LPOLY */
static double dfa;
static double plg[4][9];
static double ctloc, stloc;
static double c2tloc, s2tloc;
static double s3tloc, c3tloc;
static double apdf, apt[4];

void glatf(double lat, double &gv, double &reff) noexcept {
  constexpr double dgtr = 1.74533e-2;
  const double c2 = std::cos(2e0*dgtr*lat);
  gv = 980.616e0 * (1e0 - 0.0026373e0 * c2);
  reff = 2e0 * gv / (3.085462e-6 + 2.27e-9 * c2) * 1.0e-5;
}

double ccor(double alt, double r, double h1, double zh) noexcept {
  /// CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
  /// ALT - altitude
  /// R - target ratio
  /// H1 - transition scale length
  /// ZH - altitude of 1/2 R
  const double e = (alt - zh) / h1;
  if (e > 70)
    return std::exp(0);
  if (e < -70)
    return std::exp(r);
  
  const double ex = std::exp(e);
  // e = r / (1e0 + ex);
  return std::exp(r / (1e0+ex));
}

double ccor2(double alt, double r, double h1, double zh, double h2) noexcept {
  /// CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
  /// ALT - altitude
  /// R - target ratio
  /// H1 - transition scale length
  /// ZH - altitude of 1/2 R
  /// H2 - transition scale length #2 ?
  const double e1 = (alt - zh) / h1;
  const double e2 = (alt - zh) / h2;
  if ((e1 > 70) || (e2 > 70))
    return std::exp(0);
  if ((e1 < -70) && (e2 < -70))
    return std::exp(r);
  const double ex1 = std::exp(e1);
  const double ex2 = std::exp(e2);
  const double ccor2v = r / (1e0 + 0.5e0 * (ex1 + ex2));
  return std::exp(ccor2v);
}

double scalh(double alt, double xm, double temp) noexcept {
  constexpr double rgas = 831.4e0;
  double g = gsurf / (std::pow((1e0 + alt / re), 2e0));
  g = rgas * temp / (g * xm);
  return g;
}

double dnet(double dd, double dm, double zhm, double xmm, double xm) noexcept {
  /// TURBOPAUSE CORRECTION FOR MSIS MODELS
  ///    Root mean density
  ///     DD - diffusive density
  ///     DM - full mixed density
  ///     ZHM - transition scale length
  ///     XMM - full mixed molecular weight
  ///     XM  - species molecular weight
  ///     DNET - combined density
  double a = zhm / (xmm - xm);
  if (!((dm > 0) && (dd > 0))) {
    fprintf(stderr, "[ERROR] dnet log error %e %e %e (traceback: %s)\n", dm, dd,
            xm, __func__);
    if ((dd == 0) && (dm == 0))
      dd = 1;
    if (dm == 0)
      return dd;
    if (dd == 0)
      return dm;
  }
  double ylog = a * log(dm / dd);
  if (ylog < -10)
    return dd;
  if (ylog > 10)
    return dm;
  a = dd * std::pow((1e0 + std::exp(ylog)), (1e0 / a));
  return a;
}

void splini(double *xa, double *ya, double *y2a, int n, double x, double &y) noexcept {
  /// INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
  ///    XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
  ///    Y2A: ARRAY OF SECOND DERIVATIVES
  ///    N: SIZE OF ARRAYS XA,YA,Y2A
  ///    X: ABSCISSA ENDPOINT FOR INTEGRATION
  ///    Y: OUTPUT VALUE
  double yi = 0;
  int klo = 0;
  int khi = 1;
  double xx, h, a, b, a2, b2;
  while ((x > xa[klo]) && (khi < n)) {
    xx = x;
    if (khi < (n - 1)) {
      if (x < xa[khi])
        xx = x;
      else
        xx = xa[khi];
    }
    h = xa[khi] - xa[klo];
    a = (xa[khi] - xx) / h;
    b = (xx - xa[klo]) / h;
    a2 = a * a;
    b2 = b * b;
    yi += ((1e0 - a2) * ya[klo] / 2e0 + b2 * ya[khi] / 2e0 +
           ((-(1e0 + a2 * a2) / 4e0 + a2 / 2e0) * y2a[klo] +
            (b2 * b2 / 4e0 - b2 / 2e0) * y2a[khi]) *
               h * h / 6e0) *
          h;
    klo++;
    khi++;
  }
  y = yi;
}

void splint(double *xa, double *ya, double *y2a, int n, double x, double &y) noexcept {
  /// CALCULATE CUBIC SPLINE INTERP VALUE
  ///  ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL.
  ///  XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
  ///  Y2A: ARRAY OF SECOND DERIVATIVES
  ///  N: SIZE OF ARRAYS XA,YA,Y2A
  ///  X: ABSCISSA FOR INTERPOLATION
  ///  Y: OUTPUT VALUE
  int klo = 0;
  int khi = n - 1;
  int k;
  double h, a, b, yi;
  while ((khi - klo) > 1) {
    k = (khi + klo) / 2;
    if (xa[k] > x)
      khi = k;
    else
      klo = k;
  }
  h = xa[khi] - xa[klo];
  if (h == 0e0)
    fprintf(stderr, "[ERROR] bad XA input to splint (traceback: %s)\n", __func__);
  a = (xa[khi] - x) / h;
  b = (x - xa[klo]) / h;
  yi = a * ya[klo] + b * ya[khi] +
       ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * h * h / 6.0;
  y = yi;
}

// @todo: avoid allocation here; pass in the workspace
void spline(double *x, double *y, int n, double yp1, double ypn,
            double *y2) noexcept {
  /// CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
  /// ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL
  /// X,Y: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
  /// N: SIZE OF ARRAYS X,Y
  /// YP1,YPN: SPECIFIED DERIVATIVES AT X[0] AND X[N-1]; VALUES
  ///          >= 1E30 SIGNAL SIGNAL SECOND DERIVATIVE ZERO
  /// Y2: OUTPUT ARRAY OF SECOND DERIVATIVES
  double *u = new double[n];
  double sig, p, qn, un;
  int i, k;
  // u = malloc(sizeof(double) * (unsigned int)n);
  if (yp1 > 0.99e30) {
    y2[0] = 0;
    u[0] = 0;
  } else {
    y2[0] = -0.5e0;
    u[0] = (3e0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  }
  for (i = 1; i < (n - 1); i++) {
    sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
    p = sig * y2[i - 1] + 2e0;
    y2[i] = (sig - 1e0) / p;
    u[i] = (6e0 *
                ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
                 (y[i] - y[i - 1]) / (x[i] - x[i - 1])) /
                (x[i + 1] - x[i - 1]) -
            sig * u[i - 1]) /
           p;
  }
  if (ypn > 0.99e30) {
    qn = 0;
    un = 0;
  } else {
    qn = 0e5;
    un = (3e0 / (x[n - 1] - x[n - 2])) *
         (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
  }
  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1e0);
  for (k = n - 2; k >= 0; k--)
    y2[k] = y2[k] * y2[k + 1] + u[k];

  delete[] u;
}

inline double zeta(double zz, double zl) noexcept {
	return ((zz-zl)*(re+zl)/(re+zz));
}

double densm(double alt, double d0, double xm, double &tz, int mn3, double *zn3,
             double *tn3, double *tgn3, int mn2, double *zn2, double *tn2,
             double *tgn2) {
  // Calculate Temperature and Density Profiles for lower atmos.
  double xs[10], ys[10], y2out[10];
  double rgas = 831.4e0;
  //double z, z1, z2, t1, t2, zg, zgdif;
  //double yd1, yd2;
  //double x, y, yi;
  //double expl, gamm, glb;
  //double densm_tmp;
  //int k;
  double densm_tmp = d0;
  if (alt > zn2[0]) {
    if (xm == 0e0)
      return tz;
    else
      return d0;
  }

  // stratosphere/mesosphere temperature
  double z;
  if (alt > zn2[mn2 - 1])
    z = alt;
  else
    z = zn2[mn2 - 1];
  int mn = mn2;
  double z1 = zn2[0];
  double z2 = zn2[mn - 1];
  double t1 = tn2[0];
  double t2 = tn2[mn - 1];
  double zg = zeta(z, z1);
  double zgdif = zeta(z2, z1);

  // set up spline nodes
  for (int k = 0; k < mn; k++) {
    xs[k] = zeta(zn2[k], z1) / zgdif;
    ys[k] = 1e0 / tn2[k];
  }
  double yd1 = -tgn2[0] / (t1 * t1) * zgdif;
  double yd2 = -tgn2[1] / (t2 * t2) * zgdif * (std::pow(((re + z2) / (re + z1)), 2e0));

  // calculate spline coefficients
  spline(xs, ys, mn, yd1, yd2, y2out);
  double x = zg / zgdif;
  double y;
  splint(xs, ys, y2out, mn, x, y);

  // temperature at altitude
  tz = 1e0 / y;
  if (xm != 0e0) {
    // calaculate stratosphere / mesospehere density
    double glb = gsurf / (std::pow((1e0 + z1 / re), 2e0));
    double gamm = xm * glb * zgdif / rgas;

    // Integrate temperature profile
    double yi;
    splini(xs, ys, y2out, mn, x, yi);
    double expl = gamm * yi;
    if (expl > 50e0)
      expl = 50e0;

    // Density at altitude
    densm_tmp = densm_tmp * (t1 / tz) * std::exp(-expl);
  }

  if (alt > zn3[0]) {
    if (xm == 0e0)
      return tz;
    else
      return densm_tmp;
  }

  // troposhere / stratosphere temperature
  z = alt;
  mn = mn3;
  z1 = zn3[0];
  z2 = zn3[mn - 1];
  t1 = tn3[0];
  t2 = tn3[mn - 1];
  zg = zeta(z, z1);
  zgdif = zeta(z2, z1);

  // set up spline nodes
  for (int k = 0; k < mn; k++) {
    xs[k] = zeta(zn3[k], z1) / zgdif;
    ys[k] = 1e0 / tn3[k];
  }
  yd1 = -tgn3[0] / (t1 * t1) * zgdif;
  yd2 = -tgn3[1] / (t2 * t2) * zgdif * (pow(((re + z2) / (re + z1)), 2e0));

  // calculate spline coefficients 
  spline(xs, ys, mn, yd1, yd2, y2out);
  x = zg / zgdif;
  splint(xs, ys, y2out, mn, x, y);

  // temperature at altitude
  tz = 1e0 / y;
  if (xm != 0e0) {
    // calaculate tropospheric / stratosphere density
    double glb = gsurf / (std::pow((1e0 + z1 / re), 2e0));
    double gamm = xm * glb * zgdif / rgas;

    // Integrate temperature profile
    double yi;
    splini(xs, ys, y2out, mn, x, yi);
    double expl = gamm * yi;
    if (expl > 50e0)
      expl = 50e0;

    // Density at altitude
    densm_tmp = densm_tmp * (t1 / tz) * std::exp(-expl);
  }
  if (xm == 0e0)
    return tz;
  else
    return densm_tmp;
}

/// 3hr Magnetic activity functions
// Eq. A24d
inline double g0(double a, const double *p) noexcept {
  return (a - 4e0 +
          (p[25] - 1e0) * (a - 4e0 +
                           (std::exp(-sqrt(p[24] * p[24]) * (a - 4e0)) - 1e0) /
                               std::sqrt(p[24] * p[24])));
}

// Eq. A24c
inline double sumex(double ex) noexcept {
  return (1e0 + (1e0 - std::pow(ex, 19e0)) / (1e0 - ex) * std::pow(ex, 0.5e0));
}

// Eq. A24a
inline double sg0(double ex, const double *p, const double *ap) {
  return (g0(ap[1], p) + (g0(ap[2], p) * ex + g0(ap[3], p) * ex * ex +
                          g0(ap[4], p) * std::pow(ex, 3e0) +
                          (g0(ap[5], p) * std::pow(ex, 4e0) +
                           g0(ap[6], p) * std::pow(ex, 12e0)) *
                              (1e0 - std::pow(ex, 8e0)) / (1e0 - ex))) /
         sumex(ex);
}

double globe7(double *p, const nrlmsise_input &input,
              const nrlmsise_flags &flags) noexcept {
  // CALCULATE G(L) FUNCTION
  // Upper Thermosphere Parameters
  double t[15] = {0};
  constexpr double sr = 7.2722e-5;
  constexpr double dgtr = 1.74533e-2;
  constexpr double dr = 1.72142e-2;
  constexpr double hr = 0.2618e0;
  struct ap_array *ap;

  double tloc = input.lst;

  // calculate legendre polynomials
  double c = std::sin(input.g_lat * dgtr);
  double s = std::cos(input.g_lat * dgtr);
  double c2 = c * c;
  double c4 = c2 * c2;
  double s2 = s * s;

  double plg[4][9] plg[0][1] = c;
  plg[0][2] = 0.5 * (3.0 * c2 - 1.0);
  plg[0][3] = 0.5 * (5.0 * c * c2 - 3.0 * c);
  plg[0][4] = (35.0 * c4 - 30.0 * c2 + 3.0) / 8.0;
  plg[0][5] = (63.0 * c2 * c2 * c - 70.0 * c2 * c + 15.0 * c) / 8.0;
  plg[0][6] = (11.0 * c * plg[0][5] - 5.0 * plg[0][4]) / 6.0;
  // plg[0][7] = (13.0*c*plg[0][6] - 6.0*plg[0][5])/7.0;
  plg[1][1] = s;
  plg[1][2] = 3.0 * c * s;
  plg[1][3] = 1.5 * (5.0 * c2 - 1.0) * s;
  plg[1][4] = 2.5 * (7.0 * c2 * c - 3.0 * c) * s;
  plg[1][5] = 1.875 * (21.0 * c4 - 14.0 * c2 + 1.0) * s;
  plg[1][6] = (11.0 * c * plg[1][5] - 6.0 * plg[1][4]) / 5.0;
  // plg[1][7] = (13.0*c*plg[1][6]-7.0*plg[1][5])/6.0;
  // plg[1][8] = (15.0*c*plg[1][7]-8.0*plg[1][6])/7.0;
  plg[2][2] = 3.0 * s2;
  plg[2][3] = 15.0 * s2 * c;
  plg[2][4] = 7.5 * (7.0 * c2 - 1.0) * s2;
  plg[2][5] = 3.0 * c * plg[2][4] - 2.0 * plg[2][3];
  plg[2][6] = (11.0 * c * plg[2][5] - 7.0 * plg[2][4]) / 4.0;
  plg[2][7] = (13.0 * c * plg[2][6] - 8.0 * plg[2][5]) / 5.0;
  plg[3][3] = 15.0 * s2 * s;
  plg[3][4] = 105.0 * s2 * s * c;
  plg[3][5] = (9.0 * c * plg[3][4] - 7. * plg[3][3]) / 2.0;
  plg[3][6] = (11.0 * c * plg[3][5] - 8. * plg[3][4]) / 3.0;

  if (!(((flags.sw[7] == 0) && (flags.sw[8] == 0)) && (flags.sw[14] == 0))) {
    stloc = std::sin(hr * tloc);
    ctloc = std::cos(hr * tloc);
    s2tloc = std::sin(2.0 * hr * tloc);
    c2tloc = std::cos(2.0 * hr * tloc);
    s3tloc = std::sin(3.0 * hr * tloc);
    c3tloc = std::cos(3.0 * hr * tloc);
  }

  const double cd32 = std::cos(dr * (input.doy - p[31]));
  const double cd18 = std::cos(2.0 * dr * (input.doy - p[17]));
  const double cd14 = std::cos(dr * (input.doy - p[13]));
  const double cd39 = std::cos(2.0 * dr * (input.doy - p[38]));

  // F10.7 EFFECT
  const double df = input.f107 - input.f107A;
  dfa = input.f107A - 150.0;
  t[0] = p[19] * df * (1.0 + p[59] * dfa) + p[20] * df * df + p[21] * dfa +
         p[29] * pow(dfa, 2.0);
  const double f1 =
      1.0 + (p[47] * dfa + p[19] * df + p[20] * df * df) * flags.swc[1];
  const double f2 =
      1.0 + (p[49] * dfa + p[19] * df + p[20] * df * df) * flags.swc[1];

  //  time independent
  t[1] = (p[1] * plg[0][2] + p[2] * plg[0][4] + p[22] * plg[0][6]) +
         (p[14] * plg[0][2]) * dfa * flags.swc[1] + p[26] * plg[0][1];

  //  symmetrical annual
  t[2] = p[18] * cd32;

  //  symmetrical semiannual
  t[3] = (p[15] + p[16] * plg[0][2]) * cd18;

  //  asymmetrical annual
  t[4] = f1 * (p[9] * plg[0][1] + p[10] * plg[0][3]) * cd14;

  //  asymmetrical semiannual
  t[5] = p[37] * plg[0][1] * cd39;

  // diurnal
  if (flags.sw[7]) {
    const double t71 = (p[11] * plg[1][2]) * cd14 * flags.swc[5];
    const double t72 = (p[12] * plg[1][2]) * cd14 * flags.swc[5];
    t[6] =
        f2 * ((p[3] * plg[1][1] + p[4] * plg[1][3] + p[27] * plg[1][5] + t71) *
                  ctloc +
              (p[6] * plg[1][1] + p[7] * plg[1][3] + p[28] * plg[1][5] + t72) *
                  stloc);
  }

  // semidiurnal
  if (flags.sw[8]) {
    const double t81 =
        (p[23] * plg[2][3] + p[35] * plg[2][5]) * cd14 * flags.swc[5];
    const double t82 =
        (p[33] * plg[2][3] + p[36] * plg[2][5]) * cd14 * flags.swc[5];
    t[7] = f2 * ((p[5] * plg[2][2] + p[41] * plg[2][4] + t81) * c2tloc +
                 (p[8] * plg[2][2] + p[42] * plg[2][4] + t82) * s2tloc);
  }

  // terdiurnal
  if (flags.sw[14]) {
    t[13] = f2 * ((p[39] * plg[3][3] + (p[93] * plg[3][4] + p[46] * plg[3][6]) *
                                           cd14 * flags.swc[5]) *
                      s3tloc +
                  (p[40] * plg[3][3] + (p[94] * plg[3][4] + p[48] * plg[3][6]) *
                                           cd14 * flags.swc[5]) *
                      c3tloc);
  }

  // magnetic activity based on daily ap
  if (flags.sw[9] == -1) {
    ap = input.ap_a;
    if (p[51] != 0) {
      double exp1 = std::exp(
          -10800.0 * std::sqrt(p[51] * p[51]) /
          (1.0 + p[138] * (45.0 - std::sqrt(input.g_lat * input.g_lat))));
      if (exp1 > 0.99999)
        exp1 = 0.99999;
      if (p[24] < 1.0e-4)
        p[24] = 1.0e-4;
      apt[0] = sg0(exp1, p, ap.a);
      /* apt[1]=sg2(exp1,p,ap->a);
         apt[2]=sg0(exp2,p,ap->a);
         apt[3]=sg2(exp2,p,ap->a);
      */
      if (flags.sw[9]) {
        t[8] = apt[0] *
               (p[50] + p[96] * plg[0][2] + p[54] * plg[0][4] +
                (p[125] * plg[0][1] + p[126] * plg[0][3] + p[127] * plg[0][5]) *
                    cd14 * flags.swc[5] +
                (p[128] * plg[1][1] + p[129] * plg[1][3] + p[130] * plg[1][5]) *
                    flags.swc[7] * std::cos(hr * (tloc - p[131])));
      }
    }
  } else {
    double apd = input.ap - 4.0;
    double p44 = p[43];
    double p45 = p[44];
    if (p44 < 0)
      p44 = 1.0e-5;
    apdf = apd + (p45 - 1.0) * (apd + (std::exp(-p44 * apd) - 1.0) / p44);
    if (flags.sw[9]) {
      t[8] = apdf *
             (p[32] + p[45] * plg[0][2] + p[34] * plg[0][4] +
              (p[100] * plg[0][1] + p[101] * plg[0][3] + p[102] * plg[0][5]) *
                  cd14 * flags.swc[5] +
              (p[121] * plg[1][1] + p[122] * plg[1][3] + p[123] * plg[1][5]) *
                  flags.swc[7] * cos(hr * (tloc - p[124])));
    }
  }

  if ((flags.sw[10]) && (input.g_long > -1000.0)) {

    // longitudinal
    if (flags.sw[11]) {
      t[10] =
          (1.0 + p[80] * dfa * flags.swc[1]) *
          ((p[64] * plg[1][2] + p[65] * plg[1][4] + p[66] * plg[1][6] +
            p[103] * plg[1][1] + p[104] * plg[1][3] + p[105] * plg[1][5] +
            flags.swc[5] *
                (p[109] * plg[1][1] + p[110] * plg[1][3] + p[111] * plg[1][5]) *
                cd14) *
               std::cos(dgtr * input.g_long) +
           (p[90] * plg[1][2] + p[91] * plg[1][4] + p[92] * plg[1][6] +
            p[106] * plg[1][1] + p[107] * plg[1][3] + p[108] * plg[1][5] +
            flags.swc[5] *
                (p[112] * plg[1][1] + p[113] * plg[1][3] + p[114] * plg[1][5]) *
                cd14) *
               std::sin(dgtr * input.g_long));
    }

    // ut and mixed ut, longitude
    if (flags.sw[12]) {
      t[11] = (1.0 + p[95] * plg[0][1]) * (1.0 + p[81] * dfa * flags.swc[1]) *
              (1.0 + p[119] * plg[0][1] * flags.swc[5] * cd14) *
              ((p[68] * plg[0][1] + p[69] * plg[0][3] + p[70] * plg[0][5]) *
               std::cos(sr * (input.sec - p[71])));
      t[11] += flags.swc[11] *
               (p[76] * plg[2][3] + p[77] * plg[2][5] + p[78] * plg[2][7]) *
               std::cos(sr * (input.sec - p[79]) + 2.0 * dgtr * input.g_long) *
               (1.0 + p[137] * dfa * flags.swc[1]);
    }

    // ut, longitude magnetic activity
    if (flags.sw[13]) {
      if (flags.sw[9] == -1) {
        if (p[51]) {
          t[12] =
              apt[0] * flags.swc[11] * (1. + p[132] * plg[0][1]) *
                  ((p[52] * plg[1][2] + p[98] * plg[1][4] + p[67] * plg[1][6]) *
                   std::cos(dgtr * (input.g_long - p[97]))) +
              apt[0] * flags.swc[11] * flags.swc[5] *
                  (p[133] * plg[1][1] + p[134] * plg[1][3] +
                   p[135] * plg[1][5]) *
                  cd14 * std::cos(dgtr * (input.g_long - p[136])) +
              apt[0] * flags.swc[12] *
                  (p[55] * plg[0][1] + p[56] * plg[0][3] + p[57] * plg[0][5]) *
                  std::cos(sr * (input.sec - p[58]));
        }
      } else {
        t[12] =
            apdf * flags.swc[11] * (1.0 + p[120] * plg[0][1]) *
                ((p[60] * plg[1][2] + p[61] * plg[1][4] + p[62] * plg[1][6]) *
                 std::cos(dgtr * (input.g_long - p[63]))) +
            apdf * flags.swc[11] * flags.swc[5] *
                (p[115] * plg[1][1] + p[116] * plg[1][3] + p[117] * plg[1][5]) *
                cd14 * std::cos(dgtr * (input.g_long - p[118])) +
            apdf * flags.swc[12] *
                (p[83] * plg[0][1] + p[84] * plg[0][3] + p[85] * plg[0][5]) *
                std::cos(sr * (input.sec - p[75]));
      }
    }
  }

  // parms not used: 82, 89, 99, 139-149
  double tinf = p[30];
  for (int i = 0; i < 14; i++)
    tinf += std::fabs(flags.sw[i + 1]) * t[i];
  return tinf;
}

double glob7s(double *p, const nrlmsise_input &input,
              const nrlmsise_flags &flags) noexcept {
  // version of globe for lower atmosphere 10/26/99
  const double pset = 2.0;
  double t[14] = {0};
  const double dr = 1.72142e-2;
  const double dgtr = 1.74533e-2;
  // confirm parameter set
  if (p[99] == 0)
    p[99] = pset;
  if (p[99] != pset) {
    fprintf(stderr, "[ERROR] Wrong parameter set for glob7s (traceback: %s)\n",
            __func__);
    return -1;
  }
  const double cd32 = std::cos(dr * (input.doy - p[31]));
  const double cd18 = std::cos(2.0 * dr * (input.doy - p[17]));
  const double cd14 = std::cos(dr * (input.doy - p[13]));
  const double cd39 = std::cos(2.0 * dr * (input.doy - p[38]));

  // F10.7
  t[0] = p[21] * dfa;

  // time independent
  t[1] = p[1] * plg[0][2] + p[2] * plg[0][4] + p[22] * plg[0][6] +
         p[26] * plg[0][1] + p[14] * plg[0][3] + p[59] * plg[0][5];

  // symmetrical annual
  t[2] = (p[18] + p[47] * plg[0][2] + p[29] * plg[0][4]) * cd32;

  // symmetrical semiannual
  t[3] = (p[15] + p[16] * plg[0][2] + p[30] * plg[0][4]) * cd18;

  // asymmetrical annual
  t[4] = (p[9] * plg[0][1] + p[10] * plg[0][3] + p[20] * plg[0][5]) * cd14;

  // asymmetrical semiannual
  t[5] = (p[37] * plg[0][1]) * cd39;

  // diurnal
  if (flags.sw[7]) {
    double t71 = p[11] * plg[1][2] * cd14 * flags.swc[5];
    double t72 = p[12] * plg[1][2] * cd14 * flags.swc[5];
    t[6] = ((p[3] * plg[1][1] + p[4] * plg[1][3] + t71) * ctloc +
            (p[6] * plg[1][1] + p[7] * plg[1][3] + t72) * stloc);
  }

  // semidiurnal
  if (flags.sw[8]) {
    double t81 = (p[23] * plg[2][3] + p[35] * plg[2][5]) * cd14 * flags.swc[5];
    double t82 = (p[33] * plg[2][3] + p[36] * plg[2][5]) * cd14 * flags.swc[5];
    t[7] = ((p[5] * plg[2][2] + p[41] * plg[2][4] + t81) * c2tloc +
            (p[8] * plg[2][2] + p[42] * plg[2][4] + t82) * s2tloc);
  }

  // terdiurnal
  if (flags.sw[14]) {
    t[13] = p[39] * plg[3][3] * s3tloc + p[40] * plg[3][3] * c3tloc;
  }

  // magnetic activity
  if (flags.sw[9]) {
    if (flags.sw[9] == 1)
      t[8] = apdf * (p[32] + p[45] * plg[0][2] * flags.swc[2]);
    if (flags.sw[9] == -1)
      t[8] = (p[50] * apt[0] + p[96] * plg[0][2] * apt[0] * flags.swc[2]);
  }

  // longitudinal
  if (!((flags.sw[10] == 0) || (flags.sw[11] == 0) ||
        (input.g_long <= -1000.0))) {
    t[10] =
        (1.0 +
         plg[0][1] *
             (p[80] * flags.swc[5] * std::cos(dr * (input.doy - p[81])) +
              p[85] * flags.swc[6] * std::cos(2.0 * dr * (input.doy - p[86]))) +
         p[83] * flags.swc[3] * std::cos(dr * (input.doy - p[84])) +
         p[87] * flags.swc[4] * std::cos(2.0 * dr * (input.doy - p[88]))) *
        ((p[64] * plg[1][2] + p[65] * plg[1][4] + p[66] * plg[1][6] +
          p[74] * plg[1][1] + p[75] * plg[1][3] + p[76] * plg[1][5]) *
             std::cos(dgtr * input.g_long) +
         (p[90] * plg[1][2] + p[91] * plg[1][4] + p[92] * plg[1][6] +
          p[77] * plg[1][1] + p[78] * plg[1][3] + p[79] * plg[1][5]) *
             std::sin(dgtr * input.g_long));
  }
  double tt = 0;
  for (int i = 0; i < 14; i++)
    tt += fabs(flags.sw[i + 1]) * t[i];
  return tt;
}

void gtd7(const nrlmsise_input &input, const nrlmsise_flags &flags,
          nrlmsise_output &output) noexcept {
  const int mn3 = 5;
  const double zn3[5] = {32.5, 20.0, 15.0, 10.0, 0.0};
  const int mn2 = 4;
  const double zn2[4] = {72.5, 55.0, 45.0, 32.5};
  const double zmix = 62.5;
  nrlmsise_output soutput;

  tselec(flags);

  // Latitude variation of gravity (none for sw[2]=0)
  double xlat = input.g_lat;
  if (flags.sw[2] == 0)
    xlat = 45.0;
  glatf(xlat, gsurf, re);

  xmm = pdm[2][4];

  // thermosphere / mesosphere (above zn2[0])
  double altt;
  if (input.alt > zn2[0])
    altt = input.alt;
  else
    altt = zn2[0];

  double tmp = input.alt;
  input.alt = altt;
  gts7(input, flags, &soutput);
  altt = input.alt;
  input.alt = tmp;
  double dm28m;
  if (flags.sw[0]) // metric adjustment
    dm28m = dm28 * 1.0E6;
  else
    dm28m = dm28;
  output.t[0] = soutput.t[0];
  output.t[1] = soutput.t[1];
  if (input.alt >= zn2[0]) {
    for (int i = 0; i < 9; i++)
      output.d[i] = soutput.d[i];
    return;
  }

  // Lower mesosphere/upper stratosphere (between zn3[0] and zn2[0])
  // Temperature at nodes and gradients at end nodes
  // Inverse temperature a linear function of spherical harmonics
  meso_tgn2[0] = meso_tgn1[1];
  meso_tn2[0] = meso_tn1[4];
  meso_tn2[1] = pma[0][0] * pavgm[0] /
                (1.0 - flags.sw[20] * glob7s(pma[0], input, flags));
  meso_tn2[2] = pma[1][0] * pavgm[1] /
                (1.0 - flags.sw[20] * glob7s(pma[1], input, flags));
  meso_tn2[3] =
      pma[2][0] * pavgm[2] /
      (1.0 - flags.sw[20] * flags.sw[22] * glob7s(pma[2], input, flags));
  meso_tgn2[1] =
      pavgm[8] * pma[9][0] *
      (1.0 + flags.sw[20] * flags.sw[22] * glob7s(pma[9], input, flags)) *
      meso_tn2[3] * meso_tn2[3] / (pow((pma[2][0] * pavgm[2]), 2.0));
  meso_tn3[0] = meso_tn2[3];

  if (input.alt <= zn3[0]) {
    // Lower stratosphere and troposphere (below zn3[0])
    // Temperature at nodes and gradients at end nodes
    // Inverse temperature a linear function of spherical harmonics

    meso_tgn3[0] = meso_tgn2[1];
    meso_tn3[1] = pma[3][0] * pavgm[3] /
                  (1.0 - flags.sw[22] * glob7s(pma[3], input, flags));
    meso_tn3[2] = pma[4][0] * pavgm[4] /
                  (1.0 - flags.sw[22] * glob7s(pma[4], input, flags));
    meso_tn3[3] = pma[5][0] * pavgm[5] /
                  (1.0 - flags.sw[22] * glob7s(pma[5], input, flags));
    meso_tn3[4] = pma[6][0] * pavgm[6] /
                  (1.0 - flags.sw[22] * glob7s(pma[6], input, flags));
    meso_tgn3[1] = pma[7][0] * pavgm[7] *
                   (1.0 + flags.sw[22] * glob7s(pma[7], input, flags)) *
                   meso_tn3[4] * meso_tn3[4] /
                   (pow((pma[6][0] * pavgm[6]), 2.0));
  }

  // linear transition to full mixing below zn2[0]
  double dmc = 0;
  if (input.alt > zmix)
    dmc = 1.0 - (zn2[0] - input.alt) / (zn2[0] - zmix);
  const double dz28 = soutput.d[2];

  // N2 density
  const double dmr = soutput.d[2] / dm28m - 1.0;
  double xmm, tz;
  output.d[2] = densm(input.alt, dm28m, xmm, tz, mn3, zn3, meso_tn3, meso_tgn3,
                      mn2, zn2, meso_tn2, meso_tgn2);
  output.d[2] = output.d[2] * (1.0 + dmr * dmc);

  // HE density
  dmr = soutput.d[0] / (dz28 * pdm[0][1]) - 1.0;
  output.d[0] = output.d[2] * pdm[0][1] * (1.0 + dmr * dmc);

  // O density
  output.d[1] = 0;
  output.d[8] = 0;

  // O2 density
  dmr = soutput.d[3] / (dz28 * pdm[3][1]) - 1.0;
  output.d[3] = output.d[2] * pdm[3][1] * (1.0 + dmr * dmc);

  // AR density
  dmr = soutput.d[4] / (dz28 * pdm[4][1]) - 1.0;
  output.d[4] = output.d[2] * pdm[4][1] * (1.0 + dmr * dmc);

  // Hydrogen density
  output.d[6] = 0;

  // Atomic nitrogen density
  output.d[7] = 0;

  // Total mass density
  output.d[5] =
      1.66e-24 * (4.0 * output.d[0] + 16.0 * output.d[1] + 28.0 * output.d[2] +
                  32.0 * output.d[3] + 40.0 * output.d[4] + output.d[6] +
                  14.0 * output.d[7]);

  if (flags.sw[0])
    output.d[5] = output.d[5] / 1000;

  // temperature at altitude
  dd = densm(input.alt, 1.0, 0, tz, mn3, zn3, meso_tn3, meso_tgn3, mn2, zn2,
             meso_tn2, meso_tgn2);
  output.t[1] = tz;
}

void gtd7d(const nrlmsise_input &input, const nrlmsise_flags &flags,
           nrlmsise_output &output) noexcept {
  gtd7(input, flags, output);
  output.d[5] =
      1.66e-24 * (4.0 * output.d[0] + 16.0 * output.d[1] + 28.0 * output.d[2] +
                  32.0 * output.d[3] + 40.0 * output.d[4] + output.d[6] +
                  14.0 * output.d[7] + 16.0 * output.d[8]);
  if (flags.sw[0])
    output.d[5] = output.d[5] / 1000;
}

void ghp7(const nrlmsise_input &input, const nrlmsise_flags &flags,
          nrlmsise_output &output, double press) noexcept {
  double bm = 1.3806e-19;
  double rgas = 831.4;
  double test = 0.00043;
  double ltest = 12;

  const double pl = std::log10(press);
  double zi, z;
  if (pl >= -5.0) {
    if (pl > 2.5)
      zi = 18.06 * (3.00 - pl);
    else if ((pl > 0.075) && (pl <= 2.5))
      zi = 14.98 * (3.08 - pl);
    else if ((pl > -1) && (pl <= 0.075))
      zi = 17.80 * (2.72 - pl);
    else if ((pl > -2) && (pl <= -1))
      zi = 14.28 * (3.64 - pl);
    else if ((pl > -4) && (pl <= -2))
      zi = 12.72 * (4.32 - pl);
    else
      zi = 25.3 * (0.11 - pl);
    double cl = input.g_lat / 90.0;
    double cl2 = cl * cl;
    double cd;
    if (input.doy < 182)
      cd = (1.0 - (double)input.doy) / 91.25;
    else
      cd = ((double)input.doy) / 91.25 - 3.0;
    double ca = 0;
    if ((pl > -1.11) && (pl <= -0.23))
      ca = 1.0;
    if (pl > -0.23)
      ca = (2.79 - pl) / (2.79 + 0.23);
    if ((pl <= -1.11) && (pl > -3))
      ca = (-2.93 - pl) / (-2.93 + 1.11);
    z = zi - 4.87 * cl * cd * ca - 1.64 * cl2 * ca + 0.31 * ca * cl;
  } else
    z = 22.0 * pow((pl + 4.0), 2.0) + 110.0;

  // iteration  loop
  int l = 0;
  do {
    l++;
    input.alt = z;
    gtd7(input, flags, output);
    z = input.alt;
    double xn = output.d[0] + output.d[1] + output.d[2] + output.d[3] +
                output.d[4] + output.d[6] + output.d[7];
    double p = bm * xn * output.t[1];
    if (flags.sw[0])
      p = p * 1.0E-6;
    double diff = pl - log10(p);
    if (sqrt(diff * diff) < test)
      return;
    if (l == ltest) {
      fprintf(
          stderr,
          "ERROR: ghp7 not converging for press %e, diff %e (traceback: %s)\n",
          press, diff, __func__);
      return;
    }
    double xm = output.d[5] / xn / 1.66e-24;
    if (flags.sw[0])
      xm = xm * 1.0e3;
    double g = gsurf / (std::pow((1.0 + z / re), 2.0));
    double sh = rgas * output.t[1] / (xm * g);

    // new altitude estimate using scale height
    if (l < 6)
      z = z - sh * diff * 2.302;
    else
      z = z - sh * diff;
  } while (1 == 1);
}

void gts7(const nrlmsise_input &input, struct nrlmsise_flags *flags,
          struct nrlmsise_output *output) {
  // Thermospheric portion of NRLMSISE-00
  // See GTD7 for more extensive comments
  // alt > 72.5 km!

  double zn1[5] = {120.0, 110.0, 100.0, 90.0, 72.5};
  int mn1 = 5;
  double dgtr = 1.74533E-2;
  double dr = 1.72142E-2;
  double alpha[9] = {-0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0};
  double altl[8] = {200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0};
  za = pdl[1][15];
  zn1[0] = za;
  for (int j = 0; j < 9; j++)
    output.d[j] = 0;

  // tinf variations not important below za or zn1(1)
  if (input.alt > zn1[0])
    tinf = ptm[0] * pt[0] * (1.0 + flags.sw[16] * globe7(pt, input, flags));
  else
    tinf = ptm[0] * pt[0];
  output.t[0] = tinf;

  //  gradient variations not important below zn1(5)
  if (input.alt > zn1[4])
    g0 = ptm[3] * ps[0] * (1.0 + flags.sw[19] * globe7(ps, input, flags));
  else
    g0 = ptm[3] * ps[0];
  tlb = ptm[1] * (1.0 + flags.sw[17] * globe7(pd[3], input, flags)) * pd[3][0];
  s = g0 / (tinf - tlb);

  // Lower thermosphere temp variations not significant for
  //  density above 300 km
  if (input.alt < 300.0) {
    meso_tn1[1] = ptm[6] * ptl[0][0] /
                  (1.0 - flags.sw[18] * glob7s(ptl[0], input, flags));
    meso_tn1[2] = ptm[2] * ptl[1][0] /
                  (1.0 - flags.sw[18] * glob7s(ptl[1], input, flags));
    meso_tn1[3] = ptm[7] * ptl[2][0] /
                  (1.0 - flags.sw[18] * glob7s(ptl[2], input, flags));
    meso_tn1[4] =
        ptm[4] * ptl[3][0] /
        (1.0 - flags.sw[18] * flags.sw[20] * glob7s(ptl[3], input, flags));
    meso_tgn1[1] =
        ptm[8] * pma[8][0] *
        (1.0 + flags.sw[18] * flags.sw[20] * glob7s(pma[8], input, flags)) *
        meso_tn1[4] * meso_tn1[4] / (std::pow((ptm[4] * ptl[3][0]), 2.0));
  } else {
    meso_tn1[1] = ptm[6] * ptl[0][0];
    meso_tn1[2] = ptm[2] * ptl[1][0];
    meso_tn1[3] = ptm[7] * ptl[2][0];
    meso_tn1[4] = ptm[4] * ptl[3][0];
    meso_tgn1[1] = ptm[8] * pma[8][0] * meso_tn1[4] * meso_tn1[4] /
                   (std::pow((ptm[4] * ptl[3][0]), 2.0));
  }

  // N2 variation factor at Zlb
  g28 = flags.sw[21] * globe7(pd[2], input, flags);

  // variation of turbopause height
  zhf = pdl[1][24] *
        (1.0 + flags.sw[5] * pdl[0][24] * std::sin(dgtr * input.g_lat) *
                   std::cos(dr * (input.doy - pt[13])));
  output.t[0] = tinf;
  xmm = pdm[2][4];
  z = input.alt;

  // N2 density

  // Diffusive density at Zlb
  db28 = pdm[2][0] * std::exp(g28) * pd[2][0];
  // Diffusive density at Alt
  output.d[2] = densu(z, db28, tinf, tlb, 28.0, alpha[2], &output.t[1], ptm[5],
                      s, mn1, zn1, meso_tn1, meso_tgn1);
  dd = output.d[2];
  // Turbopause
  zh28 = pdm[2][2] * zhf;
  zhm28 = pdm[2][3] * pdl[1][5];
  xmd = 28.0 - xmm;
  // Mixed density at Zlb
  b28 = densu(zh28, db28, tinf, tlb, xmd, (alpha[2] - 1.0), tz, ptm[5], s, mn1,
              zn1, meso_tn1, meso_tgn1);
  if ((flags.sw[15]) && (z <= altl[2])) {
    //  Mixed density at Alt
    dm28 = densu(z, b28, tinf, tlb, xmm, alpha[2], tz, ptm[5], s, mn1, zn1,
                 meso_tn1, meso_tgn1);
    //  Net density at Alt
    output.d[2] = dnet(output.d[2], dm28, zhm28, xmm, 28.0);
  }

  // HE density

  // Density variation factor at Zlb
  g4 = flags.sw[21] * globe7(pd[0], input, flags);
  //  Diffusive density at Zlb
  db04 = pdm[0][0] * std::exp(g4) * pd[0][0];
  //  Diffusive density at Alt
  output.d[0] = densu(z, db04, tinf, tlb, 4., alpha[0], output.t[1], ptm[5], s,
                      mn1, zn1, meso_tn1, meso_tgn1);
  dd = output.d[0];
  if ((flags.sw[15]) && (z < altl[0])) {
    //  Turbopause
    zh04 = pdm[0][2];
    //  Mixed density at Zlb
    b04 = densu(zh04, db04, tinf, tlb, 4. - xmm, alpha[0] - 1., output.t[1],
                ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
    //  Mixed density at Alt
    dm04 = densu(z, b04, tinf, tlb, xmm, 0., output.t[1], ptm[5], s, mn1, zn1,
                 meso_tn1, meso_tgn1);
    zhm04 = zhm28;
    //  Net density at Alt
    output.d[0] = dnet(output.d[0], dm04, zhm04, xmm, 4.);
    //  Correction to specified mixing ratio at ground
    rl = std::log(b28 * pdm[0][1] / b04);
    zc04 = pdm[0][4] * pdl[1][0];
    hc04 = pdm[0][5] * pdl[1][1];
    //  Net density corrected at Alt
    output.d[0] = output.d[0] * ccor(z, rl, hc04, zc04);
  }

  // O density

  //  Density variation factor at Zlb
  g16 = flags.sw[21] * globe7(pd[1], input, flags);
  //  Diffusive density at Zlb
  db16 = pdm[1][0] * std::exp(g16) * pd[1][0];
  //   Diffusive density at Alt
  output.d[1] = densu(z, db16, tinf, tlb, 16., alpha[1], output.t[1], ptm[5], s,
                      mn1, zn1, meso_tn1, meso_tgn1);
  dd = output.d[1];
  if ((flags.sw[15]) && (z <= altl[1])) {
    //   Turbopause
    zh16 = pdm[1][2];
    //  Mixed density at Zlb
    b16 = densu(zh16, db16, tinf, tlb, 16.0 - xmm, (alpha[1] - 1.0),
                output.t[1], ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
    //  Mixed density at Alt
    dm16 = densu(z, b16, tinf, tlb, xmm, 0., output.t[1], ptm[5], s, mn1, zn1,
                 meso_tn1, meso_tgn1);
    zhm16 = zhm28;
    //  Net density at Alt
    output.d[1] = dnet(output.d[1], dm16, zhm16, xmm, 16.);
    rl = pdm[1][1] * pdl[1][16] *
         (1.0 + flags.sw[1] * pdl[0][23] * (input.f107A - 150.0));
    hc16 = pdm[1][5] * pdl[1][3];
    zc16 = pdm[1][4] * pdl[1][2];
    hc216 = pdm[1][5] * pdl[1][4];
    output.d[1] = output.d[1] * ccor2(z, rl, hc16, zc16, hc216);
    //   Chemistry correction
    hcc16 = pdm[1][7] * pdl[1][13];
    zcc16 = pdm[1][6] * pdl[1][12];
    rc16 = pdm[1][3] * pdl[1][14];
    //  Net density corrected at Alt
    output.d[1] = output.d[1] * ccor(z, rc16, hcc16, zcc16);
  }

  // O2 density

  //   Density variation factor at Zlb
  g32 = flags.sw[21] * globe7(pd[4], input, flags);
  //  Diffusive density at Zlb
  db32 = pdm[3][0] * std::exp(g32) * pd[4][0];
  //   Diffusive density at Alt
  output.d[3] = densu(z, db32, tinf, tlb, 32., alpha[3], output.t[1], ptm[5], s,
                      mn1, zn1, meso_tn1, meso_tgn1);
  dd = output.d[3];
  if (flags.sw[15]) {
    if (z <= altl[3]) {
      //   Turbopause
      zh32 = pdm[3][2];
      //  Mixed density at Zlb
      b32 = densu(zh32, db32, tinf, tlb, 32. - xmm, alpha[3] - 1., output.t[1],
                  ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
      //  Mixed density at Alt
      dm32 = densu(z, b32, tinf, tlb, xmm, 0., output.t[1], ptm[5], s, mn1, zn1,
                   meso_tn1, meso_tgn1);
      zhm32 = zhm28;
      //  Net density at Alt
      output.d[3] = dnet(output.d[3], dm32, zhm32, xmm, 32.);
      //   Correction to specified mixing ratio at ground
      rl = log(b28 * pdm[3][1] / b32);
      hc32 = pdm[3][5] * pdl[1][7];
      zc32 = pdm[3][4] * pdl[1][6];
      output.d[3] = output.d[3] * ccor(z, rl, hc32, zc32);
    }
    //  Correction for general departure from diffusive equilibrium above Zlb
    hcc32 = pdm[3][7] * pdl[1][22];
    hcc232 = pdm[3][7] * pdl[0][22];
    zcc32 = pdm[3][6] * pdl[1][21];
    rc32 = pdm[3][3] * pdl[1][23] *
           (1. + flags.sw[1] * pdl[0][23] * (input.f107A - 150.));
    //  Net density corrected at Alt
    output.d[3] = output.d[3] * ccor2(z, rc32, hcc32, zcc32, hcc232);
  }

  // AR density

  //   Density variation factor at Zlb
  g40 = flags.sw[21] * globe7(pd[5], input, flags);
  //  Diffusive density at Zlb
  db40 = pdm[4][0] * std::exp(g40) * pd[5][0];
  //   Diffusive density at Alt
  output.d[4] = densu(z, db40, tinf, tlb, 40., alpha[4], output.t[1], ptm[5], s,
                      mn1, zn1, meso_tn1, meso_tgn1);
  dd = output.d[4];
  if ((flags.sw[15]) && (z <= altl[4])) {
    //   Turbopause
    zh40 = pdm[4][2];
    //  Mixed density at Zlb
    b40 = densu(zh40, db40, tinf, tlb, 40. - xmm, alpha[4] - 1., output.t[1],
                ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
    //  Mixed density at Alt
    dm40 = densu(z, b40, tinf, tlb, xmm, 0., output.t[1], ptm[5], s, mn1, zn1,
                 meso_tn1, meso_tgn1);
    zhm40 = zhm28;
    //  Net density at Alt
    output.d[4] = dnet(output.d[4], dm40, zhm40, xmm, 40.);
    //   Correction to specified mixing ratio at ground
    rl = log(b28 * pdm[4][1] / b40);
    hc40 = pdm[4][5] * pdl[1][9];
    zc40 = pdm[4][4] * pdl[1][8];
    //  Net density corrected at Alt
    output.d[4] = output.d[4] * ccor(z, rl, hc40, zc40);
  }

  // hydrogen density

  //   Density variation factor at Zlb
  g1 = flags.sw[21] * globe7(pd[6], input, flags);
  //  Diffusive density at Zlb
  db01 = pdm[5][0] * std::exp(g1) * pd[6][0];
  //   Diffusive density at Alt
  output.d[6] = densu(z, db01, tinf, tlb, 1., alpha[6], output.t[1], ptm[5], s,
                      mn1, zn1, meso_tn1, meso_tgn1);
  dd = output.d[6];
  if ((flags.sw[15]) && (z <= altl[6])) {
    //   Turbopause
    zh01 = pdm[5][2];
    //  Mixed density at Zlb
    b01 = densu(zh01, db01, tinf, tlb, 1. - xmm, alpha[6] - 1., output.t[1],
                ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
    //  Mixed density at Alt
    dm01 = densu(z, b01, tinf, tlb, xmm, 0., output.t[1], ptm[5], s, mn1, zn1,
                 meso_tn1, meso_tgn1);
    zhm01 = zhm28;
    //  Net density at Alt
    output.d[6] = dnet(output.d[6], dm01, zhm01, xmm, 1.);
    //   Correction to specified mixing ratio at ground
    rl = log(b28 * pdm[5][1] * std::sqrt(pdl[1][17] * pdl[1][17]) / b01);
    hc01 = pdm[5][5] * pdl[1][11];
    zc01 = pdm[5][4] * pdl[1][10];
    output.d[6] = output.d[6] * ccor(z, rl, hc01, zc01);
    //   Chemistry correction
    hcc01 = pdm[5][7] * pdl[1][19];
    zcc01 = pdm[5][6] * pdl[1][18];
    rc01 = pdm[5][3] * pdl[1][20];
    //  Net density corrected at Alt
    output.d[6] = output.d[6] * ccor(z, rc01, hcc01, zcc01);
  }

  // atomic nitrogen density

  //   Density variation factor at Zlb
  g14 = flags.sw[21] * globe7(pd[7], input, flags);
  //  Diffusive density at Zlb
  db14 = pdm[6][0] * exp(g14) * pd[7][0];
  //   Diffusive density at Alt
  output.d[7] = densu(z, db14, tinf, tlb, 14., alpha[7], output.t[1], ptm[5], s,
                      mn1, zn1, meso_tn1, meso_tgn1);
  dd = output.d[7];
  if ((flags.sw[15]) && (z <= altl[7])) {
    //   Turbopause
    zh14 = pdm[6][2];
    //  Mixed density at Zlb
    b14 = densu(zh14, db14, tinf, tlb, 14. - xmm, alpha[7] - 1., output.t[1],
                ptm[5], s, mn1, zn1, meso_tn1, meso_tgn1);
    //  Mixed density at Alt
    dm14 = densu(z, b14, tinf, tlb, xmm, 0., output.t[1], ptm[5], s, mn1, zn1,
                 meso_tn1, meso_tgn1);
    zhm14 = zhm28;
    //  Net density at Alt
    output.d[7] = dnet(output.d[7], dm14, zhm14, xmm, 14.);
    //   Correction to specified mixing ratio at ground
    rl = log(b28 * pdm[6][1] * std::sqrt(pdl[0][2] * pdl[0][2]) / b14);
    hc14 = pdm[6][5] * pdl[0][1];
    zc14 = pdm[6][4] * pdl[0][0];
    output.d[7] = output.d[7] * ccor(z, rl, hc14, zc14);
    //   Chemistry correction
    hcc14 = pdm[6][7] * pdl[0][4];
    zcc14 = pdm[6][6] * pdl[0][3];
    rc14 = pdm[6][3] * pdl[0][5];
    //  Net density corrected at Alt
    output.d[7] = output.d[7] * ccor(z, rc14, hcc14, zcc14);
  }

  // Anomalous oxygen density

  g16h = flags.sw[21] * globe7(pd[8], input, flags);
  db16h = pdm[7][0] * std::exp(g16h) * pd[8][0];
  tho = pdm[7][9] * pdl[0][6];
  dd = densu(z, db16h, tho, tho, 16., alpha[8], output.t[1], ptm[5], s, mn1,
             zn1, meso_tn1, meso_tgn1);
  zsht = pdm[7][5];
  zmho = pdm[7][4];
  zsho = scalh(zmho, 16.0, tho);
  output.d[8] =
      dd * std::exp(-zsht / zsho * (std::exp(-(z - zmho) / zsht) - 1.));

  // total mass density
  output.d[5] =
      1.66e-24 * (4.0 * output.d[0] + 16.0 * output.d[1] + 28.0 * output.d[2] +
                  32.0 * output.d[3] + 40.0 * output.d[4] + output.d[6] +
                  14.0 * output.d[7]);

  // temperature
  z = std::sqrt(input.alt * input.alt);
  ddum = densu(z, 1.0, tinf, tlb, 0.0, 0.0, output.t[1], ptm[5], s, mn1, zn1,
               meso_tn1, meso_tgn1);
  // (void) ddum; // silence gcc
  if (flags.sw[0]) {
    for (int i = 0; i < 9; i++)
      output.d[i] *= 1.0e6;
    output.d[5] = output.d[5] / 1000;
  }
}
