#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>

namespace spinsim {

/// Maximum number of spins. The state vector holds 2^L complex amplitudes, so
/// L = 30 already needs 16 GB. The legacy code capped this at 24 via
/// `maxLtot` in chsdpar.h; unlike the legacy code we only allocate what a run
/// actually uses, so this is a sanity bound rather than an allocation size.
inline constexpr std::size_t kMaxSpins = 30;

/// Spin index, 0-based.
///
/// Spin `s` owns bit `s` of a basis-state index: bit set means "up", bit clear
/// means "down". This matches the legacy convention documented in
/// chstepsPDDGnmr.f ("|0> <-> |down>, |1> <-> |up>") except that the legacy
/// code numbered spins from 1 and recomputed `2**(i0-1)` inside the innermost
/// loops. Config files remain 1-based; conversion happens at the io/ boundary.
using SpinIndex = std::size_t;

/// Bit mask selecting the "up" component of `spin`.
[[nodiscard]] constexpr std::uint64_t spinMask(SpinIndex spin) noexcept {
  return std::uint64_t{1} << spin;
}

/// Number of basis states for `nspins` spins.
[[nodiscard]] constexpr std::size_t stateCount(std::size_t nspins) noexcept {
  return std::size_t{1} << nspins;
}

/// An unordered pair of distinct spins, normalised so that `low < high`.
///
/// The legacy two-spin kernel carried this invariant only as a comment
/// ("c i0 > j0 !" in chstepsPDDGnmr.f:368) and relied on every call site
/// passing its arguments in the right order. Encoding it in the type means the
/// kernel cannot be called wrongly.
class SpinPair {
 public:
  constexpr SpinPair(SpinIndex a, SpinIndex b) noexcept
      : low_(a < b ? a : b), high_(a < b ? b : a) {
    assert(a != b && "a spin cannot couple to itself");
  }

  [[nodiscard]] constexpr SpinIndex low() const noexcept { return low_; }
  [[nodiscard]] constexpr SpinIndex high() const noexcept { return high_; }

  friend constexpr bool operator==(SpinPair, SpinPair) noexcept = default;

 private:
  SpinIndex low_;
  SpinIndex high_;
};

}  // namespace spinsim
