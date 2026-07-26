#pragma once

#include <array>
#include <complex>
#include <string>
#include <string_view>
#include <vector>

#include "spinsim/core/Bits.hpp"
#include "spinsim/core/StateVector.hpp"

namespace spinsim {

/// A unit vector naming a rotation axis.
struct Axis3 {
  double x = 0.0;
  double y = 0.0;
  double z = 1.0;

  /// Throws std::invalid_argument if the vector has zero length.
  [[nodiscard]] Axis3 normalized() const;
};

/// An instantaneous, ideal rotation of one or more spins:
///     U = exp(-i (theta/2) (n . sigma))
/// applied independently to each target.
///
/// Replaces the fifteen hard-coded `ispecflag` branches of
/// chspecstepsPDDGnmr.f, which spelled out each rotation as literal arithmetic
/// on state-vector indices and selected between them by magic number.
///
/// Worth knowing: that file is dead code in every run archived here --
/// `SPECIAL PROPAGATION`, the only keyword that reaches it, appears in none of
/// the .dat files. The pulses in the CPMG and Rabi sequences are real-time
/// evolution under a large transverse field (Hx = 10700 in nmrtst2eP1.dat), so
/// they have finite duration and finite error. That is the physics the work is
/// about; these ideal gates are the idealised comparison.
class SpinRotation {
 public:
  SpinRotation(Axis3 axis, double angle, std::vector<SpinIndex> targets);

  /// Rotates every spin in the system by the same angle about the same axis.
  static SpinRotation onAllSpins(Axis3 axis, double angle,
                                 std::size_t nspins);

  void apply(StateVector& psi) const;

  [[nodiscard]] const std::vector<SpinIndex>& targets() const noexcept {
    return targets_;
  }
  [[nodiscard]] double angle() const noexcept { return angle_; }
  [[nodiscard]] Axis3 axis() const noexcept { return axis_; }

 private:
  Axis3 axis_;
  double angle_;
  std::vector<SpinIndex> targets_;
};

/// An arbitrary 4x4 unitary on a pair of spins.
///
/// Covers the legacy EPR gates (ispecflag 11 through 15), which are two-spin
/// operations no single-spin rotation reproduces.
///
/// Ordering of the basis is (down-down, low-up, high-up, up-up), matching the
/// bit layout: index = (high bit)(low bit).
class TwoSpinGate {
 public:
  using Matrix = std::array<std::complex<double>, 16>;

  TwoSpinGate(Matrix matrix, SpinPair pair);

  void apply(StateVector& psi) const;

  /// Largest deviation of U^dagger U from the identity. Zero for a true
  /// unitary; used by the tests to catch normalisation slips.
  [[nodiscard]] double unitarityDefect() const;

 private:
  Matrix matrix_;
  SpinPair pair_;
};

/// Builds a rotation from a configuration string.
///
/// Accepted forms, angle in units of pi:
///   "pi:x"        pi rotation about x, on every spin
///   "pi/2:y"      half-pi rotation about y, on every spin
///   "pi:x:3"      pi rotation about x, on spin 3 only (1-based)
///
/// Throws std::invalid_argument otherwise.
[[nodiscard]] SpinRotation makeRotation(std::string_view spec,
                                        std::size_t nspins);

}  // namespace spinsim
