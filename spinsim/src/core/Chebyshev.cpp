#include "spinsim/core/Chebyshev.hpp"

#include <algorithm>
#include <cmath>
#include <format>

#include "spinsim/core/Bessel.hpp"

namespace spinsim {
namespace {

/// Length of the decaying run required before the tail is declared negligible.
/// The legacy scan inspected three consecutive coefficients; one or two can dip
/// by accident near a zero of J_n, three in a row cannot.
constexpr std::size_t kTailRun = 3;

/// Finds the first index i >= 1 such that coefficients i, i+1, i+2 are all at
/// or below `epsilon` and strictly decreasing, or npos if there is none.
///
/// This is the condition chstepsPDDGnmr.f:112 expressed directly. The legacy
/// version reached it by growing the order and recomputing every Bessel value
/// from scratch on each attempt; because the coefficients do not depend on the
/// order requested, one pass over a long enough sequence finds the same index.
std::size_t findTailStart(std::span<const double> a, double epsilon) {
  if (a.size() < kTailRun + 1) return std::string::npos;
  for (std::size_t i = 1; i + kTailRun - 1 < a.size(); ++i) {
    const double x = std::abs(a[i]);
    const double y = std::abs(a[i + 1]);
    const double z = std::abs(a[i + 2]);
    if (x <= epsilon && y <= epsilon && z <= epsilon && y < x && z < y) {
      return i;
    }
  }
  return std::string::npos;
}

}  // namespace

ChebyshevConvergenceError::ChebyshevConvergenceError(double spectralBound,
                                                     double tau,
                                                     double epsilon,
                                                     std::size_t attemptedOrder)
    : std::runtime_error(std::format(
          "Chebyshev expansion did not converge: |H| <= {:g}, tau = {:g}, "
          "alpha = {:g}, epsilon = {:g}, order tried = {}. Shorten the time "
          "step or raise epsilon.",
          spectralBound, tau, spectralBound * std::abs(tau), epsilon,
          attemptedOrder)),
      bound_(spectralBound),
      tau_(tau),
      epsilon_(epsilon),
      order_(attemptedOrder) {}

ChebyshevExpansion ChebyshevExpansion::build(TimeMode mode,
                                             double spectralBound, double tau,
                                             double epsilon,
                                             std::size_t maxOrder) {
  if (!(spectralBound >= 0.0)) {
    throw std::invalid_argument("spectral bound must be non-negative");
  }
  if (!(epsilon > 0.0)) {
    throw std::invalid_argument("epsilon must be positive");
  }

  ChebyshevExpansion out;
  const double alpha = spectralBound * std::abs(tau);

  // A step this short is the identity to within the requested accuracy. The
  // bound is left at 1 so that a Hamiltonian of vanishing norm cannot produce a
  // division by zero downstream -- the operator is never applied here anyway.
  if (alpha < 2.0 * epsilon) {
    out.coefficients_ = {1.0};
    out.bound_ = 1.0;
    return out;
  }

  // Ask for a generous sequence in one go. The expansion needs roughly alpha
  // terms plus a margin for the tail to decay through; the legacy preliminary
  // estimate was 1.1*alpha.
  const auto want = static_cast<std::size_t>(1.5 * alpha) + 64;
  const std::size_t count = std::min(want, maxOrder + kTailRun + 1);

  const BesselSequence seq = (mode == TimeMode::Real)
                                 ? besselJ(alpha, count)
                                 : besselI(alpha, count);

  const std::size_t tail = findTailStart(seq.values, epsilon);
  if (tail == std::string::npos || tail - 1 > maxOrder) {
    throw ChebyshevConvergenceError(spectralBound, tau, epsilon, count);
  }

  // Keep orders 0 .. tail-1: coefficient `tail` is already below epsilon and
  // falling, so truncating there is what bounds the error.
  out.coefficients_.assign(seq.values.begin(),
                           seq.values.begin() + static_cast<std::ptrdiff_t>(tail));
  out.bound_ = spectralBound;

  // The Bessel argument is |tau|, so a step backwards in time is taken care of
  // by J_n(-a) = (-1)^n J_n(a). Only the odd orders change sign.
  if (tau < 0.0) {
    for (std::size_t n = 1; n < out.coefficients_.size(); n += 2) {
      out.coefficients_[n] = -out.coefficients_[n];
    }
  }

  return out;
}

}  // namespace spinsim
