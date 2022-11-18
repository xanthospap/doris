int dso::pole_tide(const dso::datetime<dso::nanoseconds> &t) noexcept {
  return dso::pole_tide(t.to_fractional_years());
}

int dso::pole_tide(double tfyears, double xp, double yp) noexcept {
  // compute secular pole at t (IERS2010, Sec 7.1.4, Eq. 21 in milliarcseconds
  const double xs = 55e0 + 1677e-3 * (tfyears-2e3);
  const double ys = 3205e-1 + 3460e-3 * (tfyears-2e3);

  // assuming xp, yp in milliarcseconds, apply Eq. 25
  const double m1 = xp - xs;
  const double m2 = -(yp-ys);

  // apply Eq. 26 to find the displacements in millimeters
  const double Stheta = -9e0 * 
}
