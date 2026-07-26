#include "spinsim/ensemble/Ensemble.hpp"

#include <cmath>
#include <format>
#include <stdexcept>

#include "spinsim/core/SpinOperators.hpp"

namespace spinsim {

InitialStateKind parseInitialStateKind(std::string_view name) {
  if (name == "random") return InitialStateKind::Random;
  if (name == "up") return InitialStateKind::AllUp;
  if (name == "down") return InitialStateKind::AllDown;
  if (name == "basis") return InitialStateKind::BasisState;
  throw std::invalid_argument(std::format(
      "unknown initial state '{}'; expected random, up, down or basis", name));
}

void prepareInitialState(StateVector& psi, InitialStateKind kind, Rng& rng,
                         std::size_t basisIndex) {
  switch (kind) {
    case InitialStateKind::AllUp:
      // Every bit set: all spins up.
      psi.setBasisState(psi.size() - 1);
      return;
    case InitialStateKind::AllDown:
      psi.setBasisState(0);
      return;
    case InitialStateKind::BasisState:
      psi.setBasisState(basisIndex);
      return;
    case InitialStateKind::Random:
      break;
  }

  // A typical state: independent complex Gaussians, normalised. Its expectation
  // values approximate infinite-temperature traces to within O(2^(-L/2)), which
  // is what makes a single realisation meaningful at all.
  for (std::size_t k = 0; k < psi.size(); ++k) {
    psi.re(k) = rng.gaussian();
    psi.im(k) = rng.gaussian();
  }
  psi.normalize();
}

void prepareReferenceState(const StateVector& psi, Axis axis,
                           StateVector& reference) {
  if (reference.size() != psi.size()) reference = StateVector(psi.nspins());
  reference.setZero();
  for (SpinIndex s = 0; s < psi.nspins(); ++s) {
    switch (axis) {
      case Axis::X: axpySigmaX(1.0, psi, s, reference); break;
      case Axis::Y: axpySigmaY(1.0, psi, s, reference); break;
      case Axis::Z: axpySigmaZ(1.0, psi, s, reference); break;
    }
  }
}

std::vector<Coupling> generateDipolarLattice2D(std::size_t nspins, Rng& rng,
                                               double scale,
                                               double minSeparation) {
  if (nspins == 0) throw std::invalid_argument("need at least one spin");

  const double halfWidth = 0.5 * std::sqrt(static_cast<double>(nspins));
  std::vector<double> x(nspins, 0.0);
  std::vector<double> y(nspins, 0.0);

  std::vector<Coupling> couplings;
  couplings.reserve(nspins * (nspins - 1) / 2);

  // Spin 0 stays at the origin, as in the legacy generator.
  for (std::size_t i = 1; i < nspins; ++i) {
    double px = 0.0;
    double py = 0.0;
    bool placed = false;
    while (!placed) {
      px = halfWidth * (2.0 * rng.uniform() - 1.0);
      py = halfWidth * (2.0 * rng.uniform() - 1.0);
      placed = true;
      for (std::size_t j = 0; j < i; ++j) {
        // The legacy rejection test, kept as-is: it requires closeness in both
        // coordinates at once, so it does not enforce a minimum distance.
        if (std::abs(px - x[j]) < minSeparation &&
            std::abs(py - y[j]) < minSeparation) {
          placed = false;
          break;
        }
      }
    }
    x[i] = px;
    y[i] = py;

    for (std::size_t j = 0; j < i; ++j) {
      const double dx = px - x[j];
      const double dy = py - y[j];
      const double r2 = dx * dx + dy * dy;
      const double jz = scale / (r2 * std::sqrt(r2));  // 1/r^3
      // Secular truncation of the dipolar coupling.
      couplings.push_back(Coupling{SpinPair{j, i}, -0.5 * jz, -0.5 * jz, jz});
    }
  }
  return couplings;
}

void applyFieldDisorder(Hamiltonian& hamiltonian, Rng& rng, double sigma) {
  for (SpinIndex s = 0; s < hamiltonian.nspins(); ++s) {
    Field f = hamiltonian.field(s);
    f.z += sigma * rng.gaussian();
    hamiltonian.setField(s, f);
  }
}

Accumulator::Accumulator(std::size_t width) : mean_(width, 0.0), m2_(width, 0.0) {}

void Accumulator::add(std::span<const double> sample) {
  if (sample.size() != mean_.size()) {
    throw std::invalid_argument("sample width does not match the accumulator");
  }
  ++count_;
  const double n = static_cast<double>(count_);
  for (std::size_t k = 0; k < sample.size(); ++k) {
    const double delta = sample[k] - mean_[k];
    mean_[k] += delta / n;
    m2_[k] += delta * (sample[k] - mean_[k]);
  }
}

void Accumulator::merge(const Accumulator& other) {
  if (other.mean_.size() != mean_.size()) {
    throw std::invalid_argument("cannot merge accumulators of different width");
  }
  if (other.count_ == 0) return;
  if (count_ == 0) {
    *this = other;
    return;
  }

  const double na = static_cast<double>(count_);
  const double nb = static_cast<double>(other.count_);
  const double total = na + nb;
  for (std::size_t k = 0; k < mean_.size(); ++k) {
    const double delta = other.mean_[k] - mean_[k];
    mean_[k] += delta * nb / total;
    m2_[k] += other.m2_[k] + delta * delta * na * nb / total;
  }
  count_ += other.count_;
}

std::vector<double> Accumulator::standardDeviation() const {
  std::vector<double> out(mean_.size(), 0.0);
  if (count_ < 2) return out;
  const double denom = static_cast<double>(count_ - 1);
  for (std::size_t k = 0; k < out.size(); ++k) {
    out[k] = std::sqrt(m2_[k] / denom);
  }
  return out;
}

std::vector<double> Accumulator::standardError() const {
  std::vector<double> out = standardDeviation();
  if (count_ < 2) return out;
  const double scale = 1.0 / std::sqrt(static_cast<double>(count_));
  for (double& v : out) v *= scale;
  return out;
}

}  // namespace spinsim
