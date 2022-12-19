namespace dso {

template <typename T, int N> struct FactorialLookUpTable {
  T table[N + 1] = {(T)1, (T)1};
  constexpr FactorialLookUpTable() noexcept {
    for (int i = 2; i < N + 1; i++)
      table[i] = table[i - 1] * i;
  }

  constexpr T factorial(int i) const noexcept { return table[i]; }
  constexpr T operator()(int i) const noexcept { return table[i]; }
}; // FactorialTable

// yield (Nominator)! / (Denominator)! as T
// T should be a float, double or long double
template <typename T, int Nominator, int Denominator>
constexpr T factorialRatio() noexcept {
  constexpr int maxInt = Nominator * (Nominator >= Denominator) +
                         Denominator * (Denominator > Nominator);
  constexpr int minInt = Nominator * (Nominator <= Denominator) +
                         Denominator * (Denominator < Nominator);
  T fact = (T)1;
  for (int i = minInt + 1; i <= maxInt; i++)
    fact *= (T)i;
  return fact;
} // factorialRatio

template <typename T>
constexpr T factorialRatio(int Nominator, int Denominator) noexcept {
  const int maxInt(Nominator * (Nominator >= Denominator) +
                   Denominator * (Denominator > Nominator));
  const int minInt(Nominator * (Nominator <= Denominator) +
                   Denominator * (Denominator < Nominator));
  T fact = (T)1;
  for (int i = minInt + 1; i <= maxInt; i++)
    fact *= (T)i;
  return fact;
} // factorialRatio

} // namespace dso
