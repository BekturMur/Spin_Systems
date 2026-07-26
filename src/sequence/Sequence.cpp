#include "spinsim/sequence/Sequence.hpp"

#include <stdexcept>

namespace spinsim {
namespace {

void requireContext(const RunContext& ctx) {
  if (ctx.hamiltonian == nullptr || ctx.state == nullptr ||
      ctx.reference == nullptr || ctx.propagator == nullptr) {
    throw std::invalid_argument("run context is missing a required component");
  }
}

/// Installs the step's field table, returning the previous one so the caller
/// can put it back. An empty override leaves the system alone.
std::vector<Field> swapFields(Hamiltonian& h, const std::vector<Field>& fields) {
  if (fields.empty()) return {};
  if (fields.size() != h.nspins()) {
    throw std::invalid_argument(
        "step field table does not match the number of spins");
  }
  std::vector<Field> previous = h.fields();
  for (SpinIndex s = 0; s < h.nspins(); ++s) h.setField(s, fields[s]);
  return previous;
}

void restoreFields(Hamiltonian& h, const std::vector<Field>& previous) {
  if (previous.empty()) return;
  for (SpinIndex s = 0; s < h.nspins(); ++s) h.setField(s, previous[s]);
}

void runStep(const Step& step, RunContext& ctx);

void runEvolve(const Evolve& step, RunContext& ctx) {
  const std::vector<Field> previous = swapFields(*ctx.hamiltonian, step.fields);

  // Both the state and its companion see the same evolution; the correlator is
  // only meaningful because of that. chebsdPDDGnmr.f makes the same pair of
  // calls, once for psi and once for pz.
  ctx.propagator->evolve(*ctx.hamiltonian, step.mode, *ctx.state, step.tau,
                         step.nsteps, step.autonormalize);
  ctx.propagator->evolve(*ctx.hamiltonian, step.mode, *ctx.reference, step.tau,
                         step.nsteps, step.autonormalize);

  restoreFields(*ctx.hamiltonian, previous);
  ctx.time += step.tau * static_cast<double>(step.nsteps);
}

void runMeasure(RunContext& ctx) {
  if (ctx.recorder == nullptr) return;

  MeasureContext mc;
  mc.time = ctx.time;
  mc.state = ctx.state;
  mc.reference = ctx.reference;
  mc.hamiltonian = ctx.hamiltonian;
  mc.subsystemSize = ctx.subsystemSize;

  std::vector<double> values;
  for (const auto& observable : ctx.observables) {
    const std::size_t at = values.size();
    values.resize(at + observable->width());
    observable->measure(mc, std::span<double>(values).subspan(at));
  }
  ctx.recorder->record(ctx.time, values);
}

void runStep(const Step& step, RunContext& ctx) {
  std::visit(
      [&ctx](const auto& concrete) {
        using T = std::decay_t<decltype(concrete)>;
        if constexpr (std::is_same_v<T, Evolve>) {
          runEvolve(concrete, ctx);
        } else if constexpr (std::is_same_v<T, Pulse>) {
          concrete.rotation.apply(*ctx.state);
          concrete.rotation.apply(*ctx.reference);
        } else if constexpr (std::is_same_v<T, Measure>) {
          runMeasure(ctx);
        } else {
          for (std::size_t i = 0; i < concrete.count; ++i) {
            for (const Step& inner : concrete.body) runStep(inner, ctx);
          }
        }
      },
      step);
}

std::size_t countMeasurements(const std::vector<Step>& steps) {
  std::size_t total = 0;
  for (const Step& step : steps) {
    if (std::holds_alternative<Measure>(step)) {
      ++total;
    } else if (const auto* repeat = std::get_if<Repeat>(&step)) {
      total += repeat->count * countMeasurements(repeat->body);
    }
  }
  return total;
}

}  // namespace

void Sequence::run(RunContext& ctx) const {
  requireContext(ctx);
  for (const Step& step : steps_) runStep(step, ctx);
}

std::size_t Sequence::measurementCount() const {
  return countMeasurements(steps_);
}

std::vector<std::string> resultColumns(
    std::span<const std::unique_ptr<Observable>> observables) {
  std::vector<std::string> columns{"time"};
  for (const auto& observable : observables) {
    for (std::string& name : observable->columns()) {
      columns.push_back(std::move(name));
    }
  }
  return columns;
}

}  // namespace spinsim
