#include "spinsim/core/Hamiltonian.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "spinsim/core/SpinOperators.hpp"

namespace spinsim {
namespace {

/// Coefficients below this are treated as absent. The legacy code applied the
/// same cutoff inline (chstepsPDDGnmr.f:224) -- worth keeping, since a typical
/// run has only Hz or only Hx populated and each skipped term is a full pass
/// over 2^L amplitudes.
constexpr double kNegligible = 1.0e-14;

}  // namespace

Hamiltonian::Hamiltonian(std::size_t nspins)
    : nspins_(nspins), fields_(nspins) {
  if (nspins == 0 || nspins > kMaxSpins) {
    throw std::out_of_range("number of spins must be in [1, kMaxSpins]");
  }
}

void Hamiltonian::setField(SpinIndex spin, Field field) {
  if (spin >= nspins_) throw std::out_of_range("spin index out of range");
  fields_[spin] = field;
  diagonalValid_ = false;
}

const Field& Hamiltonian::field(SpinIndex spin) const {
  if (spin >= nspins_) throw std::out_of_range("spin index out of range");
  return fields_[spin];
}

void Hamiltonian::setCoupling(SpinPair pair, double jx, double jy, double jz) {
  if (pair.high() >= nspins_) {
    throw std::out_of_range("coupled spin index out of range");
  }
  const auto it = std::ranges::find_if(
      couplings_, [&](const Coupling& c) { return c.pair == pair; });
  if (it != couplings_.end()) {
    it->jx = jx;
    it->jy = jy;
    it->jz = jz;
  } else {
    couplings_.push_back(Coupling{pair, jx, jy, jz});
  }
  diagonalValid_ = false;
}

void Hamiltonian::rebuildDiagonal() const {
  const std::size_t nstates = stateCount(nspins_);
  diagonal_.assign(nstates, 0.0);

  // S = sigma/2 puts the field terms at hz/2 and the coupling terms at jz/4.
  // sigma_z reads +1 on a set bit and -1 on a clear one.
  for (SpinIndex s = 0; s < nspins_; ++s) {
    const double c = 0.5 * fields_[s].z;
    if (std::abs(c) <= kNegligible) continue;
    const std::size_t mask = spinMask(s);
    for (std::size_t k = 0; k < nstates; ++k) {
      diagonal_[k] += (k & mask) ? c : -c;
    }
  }
  for (const Coupling& coupling : couplings_) {
    const double c = 0.25 * coupling.jz;
    if (std::abs(c) <= kNegligible) continue;
    const std::size_t hiMask = spinMask(coupling.pair.high());
    const std::size_t loMask = spinMask(coupling.pair.low());
    for (std::size_t k = 0; k < nstates; ++k) {
      // Aligned spins give +1, anti-aligned -1.
      const bool aligned = ((k & hiMask) != 0) == ((k & loMask) != 0);
      diagonal_[k] += aligned ? c : -c;
    }
  }
  diagonalValid_ = true;
}

double Hamiltonian::spectralBound() const noexcept {
  // |S| = 1/2 bounds each field term by |H|/2, and |S x S| = 1/4 bounds each
  // coupling term by |J|/4.
  double bound = 0.0;
  for (const Field& f : fields_) {
    bound += 0.5 * (std::abs(f.x) + std::abs(f.y) + std::abs(f.z));
  }
  for (const Coupling& c : couplings_) {
    bound += 0.25 * (std::abs(c.jx) + std::abs(c.jy) + std::abs(c.jz));
  }
  return bound;
}

void Hamiltonian::accumulate(double scale, const StateVector& in,
                             StateVector& out) const {
  if (!diagonalValid_) rebuildDiagonal();

  // Everything diagonal, in one pass.
  axpyDiagonal(scale, diagonal_, in, out);

  // Transverse fields: one traversal per spin, both components at once. The
  // kernels are in Pauli matrices, so S = sigma/2 puts these at H/2.
  for (SpinIndex s = 0; s < nspins_; ++s) {
    const double cx = scale * 0.5 * fields_[s].x;
    const double cy = scale * 0.5 * fields_[s].y;
    if (std::abs(cx) > kNegligible || std::abs(cy) > kNegligible) {
      axpyTransverse(cx, cy, in, s, out);
    }
  }

  // Off-diagonal coupling channels, at J/4. sigma_x.sigma_x and sigma_y.sigma_y
  // connect |ud> <-> |du> with weight (jx + jy) and |dd> <-> |uu> with
  // (jx - jy); the latter vanishes for the secular dipolar form jx == jy, which
  // is what every 2D-NMR run here uses.
  for (const Coupling& c : couplings_) {
    const double quarter = scale * 0.25;
    const double flipFlop = quarter * (c.jx + c.jy);
    const double flipBoth = quarter * (c.jx - c.jy);
    if (std::abs(flipFlop) > kNegligible) {
      axpyFlipFlop(flipFlop, in, c.pair, out);
    }
    if (std::abs(flipBoth) > kNegligible) {
      axpyFlipBoth(flipBoth, in, c.pair, out);
    }
  }
}

void Hamiltonian::applyPhysical(const StateVector& in, StateVector& out) const {
  out.setZero();
  accumulate(1.0, in, out);
}

void Hamiltonian::applyTwiceNormalized(double bound, const StateVector& in,
                                       StateVector& out) const {
  out.setZero();
  accumulateTwiceNormalized(bound, in, out);
}

void Hamiltonian::accumulateTwiceNormalized(double bound,
                                            const StateVector& in,
                                            StateVector& out) const {
  if (!(bound > 0.0)) {
    throw std::domain_error("spectral bound must be positive");
  }
  accumulate(-2.0 / bound, in, out);
}

}  // namespace spinsim
