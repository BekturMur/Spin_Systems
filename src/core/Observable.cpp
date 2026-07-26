#include "spinsim/core/Observable.hpp"

#include <charconv>
#include <format>
#include <stdexcept>

#include "spinsim/core/SpinOperators.hpp"

namespace spinsim {
namespace {

/// out = sum_i sigma_axis(i) * in, over every spin.
void applyTotalSpin(Axis axis, const StateVector& in, StateVector& out) {
  out.setZero();
  for (SpinIndex s = 0; s < in.nspins(); ++s) {
    switch (axis) {
      case Axis::X: axpySigmaX(1.0, in, s, out); break;
      case Axis::Y: axpySigmaY(1.0, in, s, out); break;
      case Axis::Z: axpySigmaZ(1.0, in, s, out); break;
    }
  }
}

/// Grows `scratch` to match `model` if it does not already.
void fit(StateVector& scratch, const StateVector& model) {
  if (scratch.size() != model.size()) scratch = StateVector(model.nspins());
}

void requireStates(const MeasureContext& ctx, bool needsReference) {
  if (ctx.state == nullptr) {
    throw std::invalid_argument("measurement context has no state");
  }
  if (needsReference && ctx.reference == nullptr) {
    throw std::invalid_argument("this observable needs a reference state");
  }
}

}  // namespace

Axis parseAxis(std::string_view name) {
  if (name == "x" || name == "X") return Axis::X;
  if (name == "y" || name == "Y") return Axis::Y;
  if (name == "z" || name == "Z") return Axis::Z;
  throw std::invalid_argument(
      std::format("unknown axis '{}'; expected x, y or z", name));
}

std::string_view axisName(Axis axis) {
  switch (axis) {
    case Axis::X: return "x";
    case Axis::Y: return "y";
    case Axis::Z: return "z";
  }
  return "?";
}

std::string TotalSpinCorrelator::name() const {
  return std::format("corr:{}", axisName(axis_));
}

std::vector<std::string> TotalSpinCorrelator::columns() const {
  return {std::format("corr_{}_re", axisName(axis_)),
          std::format("corr_{}_im", axisName(axis_))};
}

void TotalSpinCorrelator::measure(const MeasureContext& ctx,
                                  std::span<double> out) const {
  requireStates(ctx, true);
  fit(scratch_, *ctx.state);
  applyTotalSpin(axis_, *ctx.state, scratch_);
  const auto [re, im] = scratch_.innerProduct(*ctx.reference);
  out[0] = re;
  out[1] = im;
}

std::string SingleSpinMagnetization::name() const {
  return std::format("mag:{}:{}", spin_ + 1, axisName(axis_));
}

std::vector<std::string> SingleSpinMagnetization::columns() const {
  return {std::format("mag_{}_{}", spin_ + 1, axisName(axis_))};
}

void SingleSpinMagnetization::measure(const MeasureContext& ctx,
                                      std::span<double> out) const {
  requireStates(ctx, false);
  if (spin_ >= ctx.state->nspins()) {
    throw std::out_of_range("observable names a spin outside the system");
  }
  fit(scratch_, *ctx.state);
  switch (axis_) {
    case Axis::X: applySigmaX(*ctx.state, spin_, scratch_); break;
    case Axis::Y: applySigmaY(*ctx.state, spin_, scratch_); break;
    case Axis::Z: applySigmaZ(*ctx.state, spin_, scratch_); break;
  }
  out[0] = scratch_.innerProduct(*ctx.state).first;
}

std::vector<std::string> StateNorms::columns() const {
  return {"norm_state", "norm_reference"};
}

void StateNorms::measure(const MeasureContext& ctx,
                         std::span<double> out) const {
  requireStates(ctx, true);
  out[0] = ctx.state->norm();
  out[1] = ctx.reference->norm();
}

std::vector<std::string> Energy::columns() const { return {"energy"}; }

void Energy::measure(const MeasureContext& ctx, std::span<double> out) const {
  requireStates(ctx, false);
  if (ctx.hamiltonian == nullptr) {
    throw std::invalid_argument("energy needs a Hamiltonian in the context");
  }
  fit(scratch_, *ctx.state);
  ctx.hamiltonian->applyPhysical(*ctx.state, scratch_);
  out[0] = ctx.state->innerProduct(scratch_).first;
}

namespace {

/// Splits on ':' without allocating.
std::vector<std::string_view> split(std::string_view spec) {
  std::vector<std::string_view> parts;
  std::size_t start = 0;
  while (true) {
    const std::size_t sep = spec.find(':', start);
    if (sep == std::string_view::npos) {
      parts.push_back(spec.substr(start));
      break;
    }
    parts.push_back(spec.substr(start, sep - start));
    start = sep + 1;
  }
  return parts;
}

[[noreturn]] void rejectSpec(std::string_view spec) {
  throw std::invalid_argument(std::format(
      "unknown observable '{}'; expected one of corr:<axis>, "
      "mag:<spin>:<axis>, norm, energy",
      spec));
}

}  // namespace

std::unique_ptr<Observable> makeObservable(std::string_view spec) {
  const std::vector<std::string_view> parts = split(spec);
  if (parts.empty()) rejectSpec(spec);

  if (parts[0] == "norm" && parts.size() == 1) {
    return std::make_unique<StateNorms>();
  }
  if (parts[0] == "energy" && parts.size() == 1) {
    return std::make_unique<Energy>();
  }
  if (parts[0] == "corr" && parts.size() == 2) {
    return std::make_unique<TotalSpinCorrelator>(parseAxis(parts[1]));
  }
  if (parts[0] == "mag" && parts.size() == 3) {
    // Spins are 1-based in configuration, as they were in the .dat files.
    unsigned long index = 0;
    const auto* first = parts[1].data();
    const auto* last = first + parts[1].size();
    const auto [ptr, ec] = std::from_chars(first, last, index);
    if (ec != std::errc{} || ptr != last || index == 0) {
      throw std::invalid_argument(std::format(
          "observable '{}': spin must be a positive integer", spec));
    }
    return std::make_unique<SingleSpinMagnetization>(
        static_cast<SpinIndex>(index - 1), parseAxis(parts[2]));
  }
  rejectSpec(spec);
}

}  // namespace spinsim
