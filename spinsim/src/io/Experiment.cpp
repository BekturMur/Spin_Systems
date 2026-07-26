#include "spinsim/io/Experiment.hpp"

#include <format>
#include <fstream>
#include <ostream>
#include <stdexcept>

namespace spinsim {

void runRealization(const RunConfig& config, std::size_t index,
                    VectorRecorder& out) {
  Rng rng = Rng::forRealization(config.seed, index);

  Hamiltonian hamiltonian(config.nspins);
  for (SpinIndex s = 0; s < config.nspins; ++s) {
    hamiltonian.setField(s, config.baseFields[s]);
  }

  // Order matters: every draw from `rng` must happen in the same sequence for
  // a given index, or the realisation stops being reproducible.
  if (config.generateGeometry) {
    for (const Coupling& c : generateDipolarLattice2D(config.nspins, rng,
                                                      config.dipolarScale)) {
      hamiltonian.setCoupling(c.pair, c.jx, c.jy, c.jz);
    }
  }
  if (config.fieldDisorder != 0.0) {
    applyFieldDisorder(hamiltonian, rng, config.fieldDisorder);
  }

  StateVector state(config.nspins);
  prepareInitialState(state, config.initialState, rng, config.basisIndex);

  StateVector reference(config.nspins);
  prepareReferenceState(state, config.referenceAxis, reference);

  std::vector<std::unique_ptr<Observable>> observables;
  observables.reserve(config.observables.size());
  for (const std::string& spec : config.observables) {
    observables.push_back(makeObservable(spec));
  }

  Propagator propagator(config.nspins, config.epsilon);

  RunContext ctx;
  ctx.hamiltonian = &hamiltonian;
  ctx.state = &state;
  ctx.reference = &reference;
  ctx.propagator = &propagator;
  ctx.observables = observables;
  ctx.recorder = &out;
  ctx.time = 0.0;

  config.sequence.run(ctx);
}

EnsembleResult runExperiment(const RunConfig& config, EnsembleRunner& runner) {
  std::vector<std::unique_ptr<Observable>> probes;
  probes.reserve(config.observables.size());
  for (const std::string& spec : config.observables) {
    probes.push_back(makeObservable(spec));
  }

  // resultColumns() prefixes "time", which the runner tracks separately.
  std::vector<std::string> columns = resultColumns(probes);
  columns.erase(columns.begin());

  return runner.run(config.realizations, columns,
                    [&config](std::size_t index, VectorRecorder& out) {
                      runRealization(config, index, out);
                    });
}

void writeCsv(const EnsembleResult& result, std::ostream& out) {
  const std::size_t ncols = result.columns.size();
  if (result.statistics.width() != result.times.size() * ncols) {
    throw std::runtime_error("result shape does not match its column names");
  }

  out << "time";
  for (const std::string& name : result.columns) {
    out << ',' << name << "_mean," << name << "_sd," << name << "_stderr";
  }
  out << '\n';

  const std::vector<double>& mean = result.statistics.mean();
  const std::vector<double> sd = result.statistics.standardDeviation();
  const std::vector<double> se = result.statistics.standardError();

  for (std::size_t point = 0; point < result.times.size(); ++point) {
    out << std::format("{:.10g}", result.times[point]);
    for (std::size_t c = 0; c < ncols; ++c) {
      const std::size_t k = point * ncols + c;
      out << std::format(",{:.12g},{:.12g},{:.12g}", mean[k], sd[k], se[k]);
    }
    out << '\n';
  }
}

void writeCsv(const EnsembleResult& result, const std::filesystem::path& path) {
  if (path.has_parent_path()) {
    std::filesystem::create_directories(path.parent_path());
  }
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error(
        std::format("cannot open '{}' for writing", path.string()));
  }
  writeCsv(result, out);
}

}  // namespace spinsim
