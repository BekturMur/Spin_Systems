#pragma once

#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace spinsim {

/// Whether the propagator advances real time or imaginary time.
enum class TimeMode {
  Real,       ///< exp(+i H t), expanded in J_n. Legacy `imflag = 0`.
  Imaginary,  ///< exp(H t),    expanded in I_n. Legacy `imflag = 1`.
};

/// Thrown when no expansion of acceptable length reproduces the step.
///
/// Replaces the legacy status codes -202/-203/-204, which were returned as
/// integers and had their diagnostic context smuggled out through the first
/// five slots of the results array (chstepsPDDGnmr.f:134-138).
class ChebyshevConvergenceError : public std::runtime_error {
 public:
  ChebyshevConvergenceError(double spectralBound, double tau, double epsilon,
                            std::size_t attemptedOrder);

  [[nodiscard]] double spectralBound() const noexcept { return bound_; }
  [[nodiscard]] double tau() const noexcept { return tau_; }
  [[nodiscard]] double epsilon() const noexcept { return epsilon_; }
  [[nodiscard]] std::size_t attemptedOrder() const noexcept { return order_; }

 private:
  double bound_;
  double tau_;
  double epsilon_;
  std::size_t order_;
};

/// Truncated Chebyshev expansion of the propagator over one time step.
///
/// The step operator is written as
///     sum_{n=0..order} w_n a_n T_n(Hhat),   Hhat = -H / bound,
/// where a_n are the coefficients below and w_n is a mode-dependent weight
/// applied by the propagator: w_0 = 1 and w_n = 2 (-i)^n for real time,
/// w_n = 2 (-1)^n for imaginary time.
///
/// Replaces the order-selection block duplicated between chstepsPDDGnmr.f:84-148
/// and chimstepsPDDGnmr.f, which grew the order one at a time and recomputed the
/// whole Bessel sequence on each attempt.
class ChebyshevExpansion {
 public:
  /// Builds the expansion for one step of length `tau` under a Hamiltonian
  /// whose spectrum is bounded by `spectralBound`.
  ///
  /// `epsilon` is the coefficient magnitude below which the tail is dropped
  /// (the legacy `EpsF`, defaulting to 1e-7 via chsdpar.h).
  ///
  /// Throws ChebyshevConvergenceError if the coefficients have not decayed by
  /// `maxOrder`, and std::invalid_argument for a non-positive bound or epsilon.
  static ChebyshevExpansion build(TimeMode mode, double spectralBound,
                                  double tau, double epsilon,
                                  std::size_t maxOrder = 6000);

  /// Highest retained order. Zero means the step is short enough that the
  /// identity alone represents it.
  [[nodiscard]] std::size_t order() const noexcept {
    return coefficients_.empty() ? 0 : coefficients_.size() - 1;
  }

  [[nodiscard]] std::span<const double> coefficients() const noexcept {
    return coefficients_;
  }

  /// The bound the coefficients were built against, and the one the propagator
  /// must rescale the Hamiltonian by.
  ///
  /// For a negligibly short step this is 1 rather than the true bound: the
  /// expansion is then a single constant term, the Hamiltonian is never
  /// applied, and a near-zero bound would otherwise divide by nothing. The
  /// legacy code did the same via `emax = 1.0D0 ! precaution for too small
  /// emax`.
  [[nodiscard]] double effectiveBound() const noexcept { return bound_; }

  /// True when the step was short enough to need no Hamiltonian application.
  [[nodiscard]] bool isIdentity() const noexcept { return order() == 0; }

 private:
  std::vector<double> coefficients_;
  double bound_ = 1.0;
};

}  // namespace spinsim
