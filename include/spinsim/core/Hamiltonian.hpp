#pragma once

#include <vector>

#include "spinsim/core/Bits.hpp"
#include "spinsim/core/SpinOperators.hpp"
#include "spinsim/core/StateVector.hpp"

namespace spinsim {

/// Local field acting on one spin, in the same units the legacy .dat files use.
struct Field {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

/// Anisotropic Heisenberg coupling between two spins.
struct Coupling {
  SpinPair pair;
  double jx = 0.0;
  double jy = 0.0;
  double jz = 0.0;
};

/// H = sum_i (Hx_i Sx_i + Hy_i Sy_i + Hz_i Sz_i)
///   + sum_<ij> (Jx Sx_i Sx_j + Jy Sy_i Sy_j + Jz Sz_i Sz_j),  S = sigma/2.
///
/// Replaces the eight parallel arrays the legacy code threaded through every
/// call (`Hx`, `Hy`, `Hz`, `Jx`, `Jy`, `Jz`, `asp`, `bsp`, plus `nconn` and
/// `Ltot` to say how much of each was live).
class Hamiltonian {
 public:
  explicit Hamiltonian(std::size_t nspins);

  [[nodiscard]] std::size_t nspins() const noexcept { return nspins_; }

  void setField(SpinIndex spin, Field field);
  [[nodiscard]] const Field& field(SpinIndex spin) const;

  /// Adds a coupling, or overwrites the existing one for the same pair --
  /// matching chparsfPDDG.f, where a repeated `i,j` line updates in place.
  void setCoupling(SpinPair pair, double jx, double jy, double jz);

  [[nodiscard]] const std::vector<Coupling>& couplings() const noexcept {
    return couplings_;
  }
  [[nodiscard]] const std::vector<Field>& fields() const noexcept {
    return fields_;
  }

  /// Upper bound on |H|, the sum of the operator norms of the individual terms
  /// (chstepsPDDGnmr.f:65-78). Chebyshev needs a bound, not the true extremal
  /// eigenvalue: overestimating costs expansion order, underestimating breaks
  /// convergence.
  ///
  /// Returns 0 for an empty Hamiltonian; callers propagating in time must
  /// handle that (the legacy code substituted emax = 1 in that case).
  [[nodiscard]] double spectralBound() const noexcept;

  /// out = H * in, with H the physical (spin-1/2) Hamiltonian above.
  void applyPhysical(const StateVector& in, StateVector& out) const;

  /// out = 2 * Hhat * in, where Hhat = -H / spectralBound() has spectrum
  /// inside [-1, 1].
  ///
  /// The factor of two is what the Chebyshev recurrence T_{n+1} = 2x T_n -
  /// T_{n-1} needs, so folding it in here saves a pass over the state vector
  /// per expansion order. The sign is inherited verbatim from the legacy
  /// `hefx = -Hx/emax` (chstepsPDDGnmr.f:175); Propagator's documentation
  /// records which time direction that yields, as pinned down by the test
  /// against the exact matrix exponential.
  ///
  /// `bound` must be the value returned by spectralBound(); it is passed in so
  /// the propagator computes it once per run rather than once per step.
  void applyTwiceNormalized(double bound, const StateVector& in,
                            StateVector& out) const;

  /// As applyTwiceNormalized, but adds into `out` instead of overwriting it.
  ///
  /// The Chebyshev recurrence needs exactly this: `out` arrives holding
  /// -T_{n-2} and leaves holding 2*Hhat*T_{n-1} - T_{n-2} = T_n, with no extra
  /// pass over the state vector.
  void accumulateTwiceNormalized(double bound, const StateVector& in,
                                 StateVector& out) const;

 private:
  /// out += scale * (Pauli-basis H) * in, skipping negligible coefficients.
  void accumulate(double scale, const StateVector& in, StateVector& out) const;

  /// Fills diagonal_ with the sigma_z part of the Pauli-basis operator.
  ///
  /// Every sigma_z term is diagonal in the computational basis -- the local
  /// fields and the Ising part of every coupling alike -- so all of them can be
  /// collapsed into one vector of length 2^L and applied in a single pass.
  /// Building it costs about as much as one Hamiltonian application, and it is
  /// reused across every order of the Chebyshev expansion, so a step of order
  /// sixty pays roughly two per cent for it.
  ///
  /// Rebuilt lazily: any change to a field or coupling clears the cache.
  void rebuildDiagonal() const;


  std::size_t nspins_;
  std::vector<Field> fields_;
  std::vector<Coupling> couplings_;

  mutable std::vector<double> diagonal_;
  mutable bool diagonalValid_ = false;
};

}  // namespace spinsim
