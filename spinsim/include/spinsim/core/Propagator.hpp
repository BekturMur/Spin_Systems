#pragma once

#include <cstddef>

#include "spinsim/core/Chebyshev.hpp"
#include "spinsim/core/Hamiltonian.hpp"
#include "spinsim/core/StateVector.hpp"

namespace spinsim {

/// Advances a state under a Hamiltonian by Chebyshev expansion of the step
/// operator.
///
/// Port of chstepsPDDGnmr.f (real time) and chimstepsPDDGnmr.f (imaginary
/// time), which shared everything but the Bessel family and a sign.
///
/// The Hamiltonian is supplied per call rather than held, because a pulse
/// sequence changes the fields between steps. What the propagator does own is
/// the pair of scratch state vectors the recurrence needs; the legacy routines
/// declared these as 2^24-sized locals on every call.
class Propagator {
 public:
  /// `epsilon` is the coefficient cutoff -- the legacy `EpsF`, whose default in
  /// chsdpar.h is 1e-7.
  explicit Propagator(std::size_t nspins, double epsilon = 1.0e-7);

  /// Advances `psi` by `nsteps` steps of length `tau`.
  ///
  /// The expansion is built once for the whole call, so all steps share it;
  /// this is only valid because the Hamiltonian is fixed for the duration.
  ///
  /// With `autonormalize` the state is rescaled to unit norm after each step.
  /// Real-time evolution is unitary and does not need it -- switching it on
  /// there hides any loss of accuracy rather than fixing it. Imaginary-time
  /// evolution is not norm-preserving and generally does.
  ///
  /// Throws ChebyshevConvergenceError if `tau` is too long for `epsilon`.
  void evolve(const Hamiltonian& hamiltonian, TimeMode mode, StateVector& psi,
              double tau, std::size_t nsteps, bool autonormalize = false);

  /// Expansion order used by the most recent evolve() call. Diagnostic only.
  [[nodiscard]] std::size_t lastOrder() const noexcept { return lastOrder_; }

  [[nodiscard]] double epsilon() const noexcept { return epsilon_; }

 private:
  void applyOneStep(const Hamiltonian& hamiltonian,
                    const ChebyshevExpansion& expansion, TimeMode mode,
                    StateVector& psi);

  double epsilon_;
  std::size_t lastOrder_ = 0;
  StateVector current_;   ///< T_n applied to the state at the step's start
  StateVector previous_;  ///< T_{n-1}, then reused as the recurrence target
};

}  // namespace spinsim
