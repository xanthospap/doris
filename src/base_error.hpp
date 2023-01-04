#ifndef __DORIS_BASE_STATUS_ERROR_HPP__
#define __DORIS_BASE_STATUS_ERROR_HPP__

namespace dso {
struct [[nodiscard]] iStatus {
  int status;
  constexpr iStatus(int i) noexcept : status(i){};
  explicit constexpr operator int() const noexcept { return status; }
  static constexpr iStatus ok() noexcept { return iStatus(0); }
  constexpr bool operator==(const iStatus &other) const noexcept {
    return status == other.status;
  }
  constexpr bool operator!=(const iStatus &other) const noexcept {
    return !(*this == other);
  }
};

} // namespace dso

#endif
