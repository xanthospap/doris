#ifndef __DSO_ENUM_CLASS_BITSET_HPP__
#define __DSO_ENUM_CLASS_BITSET_HPP__

#include <type_traits>
#include <vector>

namespace dso {
/// Wrapper class around an enum to make it act like a flag class.
/// The class only provides basic functionality like, construction of flags,
/// set, check and clear flags.
///
/// @tparam FlagEnum A (strongly typed) enumeration class
template <typename FlagEnum> class EnumBitset {
public:
  /// ft is the underlying enum class type (normaly int or char).
  using ft = typename std::underlying_type<FlagEnum>::type;

private:
  ft _f; ///< the flag as underlying type (ft).

public:
  /// Default constructor (no FlagEnum is set).
  EnumBitset() noexcept : _f{0} {}

  /// Constructor from a EnumBitset enumeration type (i.e. from a enumeration of
  /// the FlagEnum class).
  /// @param[in] f  A FlagEnum instance.
  explicit EnumBitset(FlagEnum f) noexcept : _f{static_cast<ft>(f)} {}

  /// Constructor from multiple EnumBitsets (via an std::initializer_list).
  /// @param[in] fs  An std::initializer_list<FlagEnum> instance; all elements
  ///                of this list will be set on on the produced instance.
  explicit EnumBitset(std::initializer_list<FlagEnum> &&fs) : _f{} {
    for (auto &f : fs) {
      _f |= static_cast<ft>(f);
    }
  }

  /// Set a FlagEnum on.
  /// @param[in] f  A FlagEnum instance; this will be "switched on".
  EnumBitset &set(FlagEnum f) noexcept {
    _f |= static_cast<ft>(f);
    return *this;
  }

  /// Clear a FlagEnum.
  /// @param[in] f  Clear (i.e. "switch off") this FlagEnum.
  EnumBitset &clear(FlagEnum f) noexcept {
    _f &= ~(static_cast<ft>(f));
    return *this;
  }

  /// Clear the instance from all FlagEnums. All FlagEnum are "switched off".
  void clear() noexcept { _f = static_cast<ft>(0); }

  /// Check if a FlagEnum is set (i.e. on).
  /// @param[in] f  Check if this FlagEnum is "switched on".
  bool check(FlagEnum f) const noexcept { return _f & static_cast<ft>(f); }

  /// Check if a EnumBitset is clean (nothing is set).
  bool is_clean() const noexcept { return !static_cast<ft>(_f); }

  /// Equality operator
  bool operator==(EnumBitset f) const noexcept { return _f == f._f; }

  /// InEquality operator
  bool operator!=(EnumBitset f) const noexcept {
    return !(this->operator==(f));
  }

}; // EnumBitset<FlagEnum>

} // namespace dso

#endif
