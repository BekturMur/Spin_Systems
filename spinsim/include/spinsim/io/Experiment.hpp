#pragma once

#include <filesystem>
#include <iosfwd>

#include "spinsim/ensemble/EnsembleRunner.hpp"
#include "spinsim/io/Config.hpp"

namespace spinsim {

/// Runs one realisation of a configured experiment.
///
/// Everything the realisation needs is derived from `index`: the spin
/// positions, the field disorder, and the initial state all come from
/// Rng::forRealization(config.seed, index). Nothing else varies between
/// realisations, so the result is a pure function of (config, index).
void runRealization(const RunConfig& config, std::size_t index,
                    VectorRecorder& out);

/// Averages the experiment over its configured ensemble.
[[nodiscard]] EnsembleResult runExperiment(const RunConfig& config,
                                           EnsembleRunner& runner);

/// Writes results as CSV: time, then mean/sd/stderr for every column.
///
/// The header names each column, taken from the observables themselves. The
/// legacy output was a bare grid of numbers whose meaning lived in the shape of
/// `vecResOut` and a parallel array of record lengths.
void writeCsv(const EnsembleResult& result, std::ostream& out);

void writeCsv(const EnsembleResult& result, const std::filesystem::path& path);

}  // namespace spinsim
