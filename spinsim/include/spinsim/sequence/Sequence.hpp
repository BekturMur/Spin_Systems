#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <variant>
#include <vector>

#include "spinsim/core/Chebyshev.hpp"
#include "spinsim/core/Gate.hpp"
#include "spinsim/core/Hamiltonian.hpp"
#include "spinsim/core/Observable.hpp"
#include "spinsim/core/Propagator.hpp"
#include "spinsim/core/StateVector.hpp"

namespace spinsim {

/// Advance time under the Hamiltonian.
///
/// `fields` replaces the per-spin fields for the duration of the step; leaving
/// it empty keeps whatever the base system carries. Couplings always come from
/// the base system, because in the 2D-NMR runs they are generated per
/// realisation rather than configured.
///
/// This is how the real pulses are expressed: a step with a large transverse
/// field and a short duration, exactly as nmrtst2eP1.dat does with Hx = 10700.
struct Evolve {
  double tau = 0.0;
  std::size_t nsteps = 1;
  TimeMode mode = TimeMode::Real;
  bool autonormalize = false;
  std::vector<Field> fields;
};

/// Apply an idealised instantaneous rotation.
struct Pulse {
  SpinRotation rotation;
};

/// Record every configured observable at the current time.
struct Measure {};

struct Repeat;

using Step = std::variant<Evolve, Pulse, Measure, Repeat>;

/// Run the body `count` times.
///
/// Nesting works, which the legacy format rejected outright: chparsgenPDDG.f
/// returns status -987 for "nested cycle", because @CYCLE was implemented as a
/// single pair of jump targets (`imark`) rather than as structure.
struct Repeat {
  std::size_t count = 1;
  std::vector<Step> body;
};

/// Where measurement results go.
class Recorder {
 public:
  virtual ~Recorder() = default;
  /// `values` is the concatenation of every observable's output, in order.
  virtual void record(double time, std::span<const double> values) = 0;
};

/// Everything one pass through a sequence needs.
///
/// Both states are propagated: `state` is |psi(t)> and `reference` is
/// A|psi(0)> carried forward through the same evolution, so that observables
/// can form the two-time correlator the echo signal is built from.
struct RunContext {
  Hamiltonian* hamiltonian = nullptr;
  StateVector* state = nullptr;
  StateVector* reference = nullptr;
  Propagator* propagator = nullptr;
  std::span<const std::unique_ptr<Observable>> observables;
  Recorder* recorder = nullptr;
  std::size_t subsystemSize = 0;
  double time = 0.0;
};

/// A pulse sequence: an ordered program of steps.
class Sequence {
 public:
  explicit Sequence(std::vector<Step> steps) : steps_(std::move(steps)) {}

  [[nodiscard]] const std::vector<Step>& steps() const noexcept {
    return steps_;
  }

  /// Executes the program, advancing `ctx.time` and both state vectors.
  void run(RunContext& ctx) const;

  /// Total number of measurement points a run will produce, counting
  /// repetitions. Lets a caller size its output buffer up front.
  [[nodiscard]] std::size_t measurementCount() const;

 private:
  std::vector<Step> steps_;
};

/// Concatenated column names for a set of observables, prefixed by "time".
[[nodiscard]] std::vector<std::string> resultColumns(
    std::span<const std::unique_ptr<Observable>> observables);

}  // namespace spinsim
