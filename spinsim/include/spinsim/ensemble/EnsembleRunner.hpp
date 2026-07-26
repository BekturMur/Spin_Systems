#pragma once

#include <cstddef>
#include <functional>
#include <span>
#include <string>
#include <vector>

#include "spinsim/ensemble/Ensemble.hpp"
#include "spinsim/sequence/Sequence.hpp"

namespace spinsim {

/// Collects measurement rows into memory, in the order the sequence produces
/// them.
class VectorRecorder final : public Recorder {
 public:
  void record(double time, std::span<const double> values) override;

  [[nodiscard]] const std::vector<double>& times() const noexcept {
    return times_;
  }
  /// Row-major: row `i` occupies [i*width, (i+1)*width).
  [[nodiscard]] const std::vector<double>& values() const noexcept {
    return values_;
  }
  [[nodiscard]] std::size_t rows() const noexcept { return times_.size(); }

  void clear();

 private:
  std::vector<double> times_;
  std::vector<double> values_;
};

/// The averaged outcome of an ensemble.
struct EnsembleResult {
  std::vector<double> times;         ///< one per measurement point
  std::vector<std::string> columns;  ///< names of the per-point values
  Accumulator statistics{0};         ///< width = times.size() * columns.size()
};

/// Runs one realisation and reports its measurements.
using RealizationJob = std::function<void(std::size_t index, VectorRecorder&)>;

/// Averages a job over many independent realisations.
///
/// The interface exists so that a cluster backend can be added without the
/// physics code knowing. The legacy program wove MPI calls directly through its
/// main routine, which is why the 2D-NMR modification had to comment out half
/// the broadcasts (chebNMR2D-v6.f:211-223) once the couplings became
/// worker-local.
class EnsembleRunner {
 public:
  virtual ~EnsembleRunner() = default;

  /// `columns` names the per-point values a job records.
  virtual EnsembleResult run(std::size_t realizations,
                             std::span<const std::string> columns,
                             const RealizationJob& job) = 0;
};

/// Runs realisations across a pool of threads in this process.
///
/// Results are bitwise independent of `threads`: jobs run in whatever order the
/// scheduler picks, but their outputs are folded into the statistics strictly
/// in realisation order. Floating-point summation is not associative, so
/// without that discipline the answer would drift with the thread count.
class LocalRunner final : public EnsembleRunner {
 public:
  /// `threads` of 0 means hardware concurrency.
  explicit LocalRunner(std::size_t threads = 0);

  EnsembleResult run(std::size_t realizations,
                     std::span<const std::string> columns,
                     const RealizationJob& job) override;

  [[nodiscard]] std::size_t threads() const noexcept { return threads_; }

 private:
  std::size_t threads_;
};

}  // namespace spinsim
