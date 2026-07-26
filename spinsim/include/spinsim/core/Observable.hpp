#pragma once

#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "spinsim/core/Bits.hpp"
#include "spinsim/core/Hamiltonian.hpp"
#include "spinsim/core/StateVector.hpp"

namespace spinsim {

/// A Cartesian axis, used to name spin components in configuration.
enum class Axis { X, Y, Z };

[[nodiscard]] Axis parseAxis(std::string_view name);
[[nodiscard]] std::string_view axisName(Axis axis);

/// Everything an observable may look at when the sequence reaches a
/// measurement point.
///
/// `state` and `reference` are both propagated through the sequence. The
/// reference starts life as A|psi_0> for some operator A, so that
/// <reference(t)| B |state(t)> is the two-time correlator
/// <psi_0| A^dagger U^dagger(t) B U(t) |psi_0>. Evaluated on a random typical
/// state this estimates the infinite-temperature correlator, which is the NMR
/// echo signal the whole program exists to compute.
struct MeasureContext {
  double time = 0.0;
  const StateVector* state = nullptr;
  const StateVector* reference = nullptr;
  const Hamiltonian* hamiltonian = nullptr;
  std::size_t subsystemSize = 0;  ///< legacy `Lsys`
};

/// A quantity recorded at each measurement point.
///
/// This interface is what removes the need to fork the source per experiment.
/// In the legacy code the choice of correlator lived in chebsdPDDGnmr.f -- line
/// 110 reads `addHSx` in CPMGZ and `addHSz` in Rabi_10_SzSz -- so measuring a
/// different component meant editing the file and copying the whole directory.
/// Roughly forty of the forty-seven directories differ in no other way.
class Observable {
 public:
  virtual ~Observable() = default;

  /// Identifier as it appears in configuration, e.g. "corr:z".
  [[nodiscard]] virtual std::string name() const = 0;

  /// Column headings, one per value emitted. Output files are self-describing
  /// as a result, unlike the legacy flat `vecResOut` array whose layout was
  /// recorded only in a side array of record lengths.
  [[nodiscard]] virtual std::vector<std::string> columns() const = 0;

  /// Writes exactly columns().size() values into `out`.
  virtual void measure(const MeasureContext& ctx,
                       std::span<double> out) const = 0;

  /// Number of values emitted.
  [[nodiscard]] std::size_t width() const { return columns().size(); }
};

/// <reference(t)| sum_i sigma_axis(i) |state(t)>, real and imaginary parts.
///
/// Together with the axis chosen when the reference state is prepared, one
/// class covers the SxSx, SySy, SzSz and SxSz variants that the legacy code
/// spread across separate directories.
class TotalSpinCorrelator final : public Observable {
 public:
  explicit TotalSpinCorrelator(Axis axis) : axis_(axis) {}

  [[nodiscard]] std::string name() const override;
  [[nodiscard]] std::vector<std::string> columns() const override;
  void measure(const MeasureContext& ctx, std::span<double> out) const override;

 private:
  Axis axis_;
  mutable StateVector scratch_;
};

/// <state(t)| sigma_axis(spin) |state(t)>, the magnetisation of one spin.
class SingleSpinMagnetization final : public Observable {
 public:
  SingleSpinMagnetization(SpinIndex spin, Axis axis)
      : spin_(spin), axis_(axis) {}

  [[nodiscard]] std::string name() const override;
  [[nodiscard]] std::vector<std::string> columns() const override;
  void measure(const MeasureContext& ctx, std::span<double> out) const override;

 private:
  SpinIndex spin_;
  Axis axis_;
  mutable StateVector scratch_;
};

/// Norms of the propagated state and of the reference state.
///
/// Under exact unitary evolution both stay at their initial values, so a drift
/// here is a direct readout of accumulated expansion error. The legacy code
/// recorded the same two numbers (chebsdPDDGnmr.f:273).
class StateNorms final : public Observable {
 public:
  [[nodiscard]] std::string name() const override { return "norm"; }
  [[nodiscard]] std::vector<std::string> columns() const override;
  void measure(const MeasureContext& ctx, std::span<double> out) const override;
};

/// Expectation of the physical Hamiltonian, <state|H|state>.
///
/// Conserved under real-time evolution; a second, independent check on the
/// propagator that the legacy code never recorded.
class Energy final : public Observable {
 public:
  [[nodiscard]] std::string name() const override { return "energy"; }
  [[nodiscard]] std::vector<std::string> columns() const override;
  void measure(const MeasureContext& ctx, std::span<double> out) const override;

 private:
  mutable StateVector scratch_;
};

/// Builds an observable from its configuration name.
///
/// Accepted forms:
///   "corr:<axis>"        total-spin correlator against the reference state
///   "mag:<spin>:<axis>"  magnetisation of one spin, `spin` being 1-based
///   "norm"               norms of the state and the reference
///   "energy"             <H>
///
/// Throws std::invalid_argument for anything else, naming what was accepted.
[[nodiscard]] std::unique_ptr<Observable> makeObservable(std::string_view spec);

}  // namespace spinsim
