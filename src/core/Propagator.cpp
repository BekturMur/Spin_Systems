#include "spinsim/core/Propagator.hpp"

#include <stdexcept>

namespace spinsim {
namespace {

/// The weight w_n multiplying a_n T_n in the expansion, excluding the factor of
/// two that every term past the zeroth carries.
///
/// Real time accumulates powers of -i; imaginary time accumulates powers of -1.
/// The legacy code tracked both with a pair of integers advanced by hand
/// (`nR`/`nI`), differently in each routine.
struct TermWeight {
  double re;
  double im;

  static TermWeight initial(TimeMode mode) {
    // The value for n = 1, since the zeroth term is handled separately.
    return mode == TimeMode::Real ? TermWeight{0.0, -1.0}
                                  : TermWeight{-1.0, 0.0};
  }

  void advance(TimeMode mode) {
    if (mode == TimeMode::Real) {
      // Multiply by -i: (a + bi)(-i) = b - ai.
      const double previousRe = re;
      re = im;
      im = -previousRe;
    } else {
      re = -re;
      im = -im;
    }
  }
};

}  // namespace

Propagator::Propagator(std::size_t nspins, double epsilon)
    : epsilon_(epsilon), current_(nspins), previous_(nspins) {
  if (!(epsilon > 0.0)) {
    throw std::invalid_argument("epsilon must be positive");
  }
}

void Propagator::applyOneStep(const Hamiltonian& hamiltonian,
                              const ChebyshevExpansion& expansion,
                              TimeMode mode, StateVector& psi) {
  const auto a = expansion.coefficients();
  const double bound = expansion.effectiveBound();

  // T_0 is the state itself; psi is then reused as the running sum.
  current_ = psi;
  psi.scale(a[0]);
  previous_.setZero();

  TermWeight w = TermWeight::initial(mode);
  for (std::size_t n = 1; n <= expansion.order(); ++n) {
    // previous_ holds T_{n-2}; negating it and adding 2*Hhat*T_{n-1} in place
    // completes the recurrence T_n = 2 Hhat T_{n-1} - T_{n-2}.
    previous_.negate();
    hamiltonian.accumulateTwiceNormalized(bound, current_, previous_);
    // T_1 = Hhat T_0, without the doubling the general recurrence carries.
    if (n == 1) previous_.scale(0.5);

    // previous_ now holds T_n and current_ holds T_{n-1}; swapping makes
    // current_ the newest term and frees the older buffer for the next round.
    previous_.swap(current_);

    psi.addScaled(2.0 * a[n] * w.re, 2.0 * a[n] * w.im, current_);
    w.advance(mode);
  }
}

void Propagator::evolve(const Hamiltonian& hamiltonian, TimeMode mode,
                        StateVector& psi, double tau, std::size_t nsteps,
                        bool autonormalize) {
  if (psi.nspins() != hamiltonian.nspins()) {
    throw std::invalid_argument("state and Hamiltonian disagree on spin count");
  }
  if (psi.nspins() != current_.nspins()) {
    throw std::invalid_argument("state does not match the propagator's size");
  }
  if (nsteps == 0) return;

  const ChebyshevExpansion expansion =
      ChebyshevExpansion::build(mode, hamiltonian.spectralBound(), tau,
                                epsilon_);
  lastOrder_ = expansion.order();

  // A step below the accuracy threshold leaves the state alone, and the
  // Hamiltonian is never applied -- which is what lets effectiveBound() be a
  // placeholder in that case.
  if (expansion.isIdentity()) return;

  for (std::size_t step = 0; step < nsteps; ++step) {
    applyOneStep(hamiltonian, expansion, mode, psi);
    if (autonormalize) psi.normalize();
  }
}

}  // namespace spinsim
