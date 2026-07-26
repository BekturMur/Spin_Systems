#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <sstream>

#include "spinsim/io/Config.hpp"
#include "spinsim/io/Experiment.hpp"

using namespace spinsim;

namespace {

/// A minimal but complete configuration, used as the base for variations.
constexpr std::string_view kMinimal = R"(
[system]
spins = 3

[field]
hz = 1.0

[measurement]
observables = ["corr:z"]

[[step]]
kind = "measure"

[[step]]
kind = "evolve"
tau = 0.1

[[step]]
kind = "measure"
)";

RunConfig parse(std::string_view text) { return parseConfig(text, "<test>"); }

}  // namespace

TEST_CASE("a minimal configuration loads", "[config]") {
  const RunConfig cfg = parse(kMinimal);
  REQUIRE(cfg.nspins == 3);
  REQUIRE(cfg.baseFields.size() == 3);
  REQUIRE(cfg.baseFields[0].z == 1.0);
  REQUIRE(cfg.sequence.measurementCount() == 2);
  REQUIRE(cfg.initialState == InitialStateKind::Random);
}

TEST_CASE("a scalar field applies to every spin", "[config]") {
  const RunConfig cfg = parse(kMinimal);
  for (const Field& f : cfg.baseFields) REQUIRE(f.z == 1.0);
}

TEST_CASE("a per-spin field array is honoured", "[config]") {
  const RunConfig cfg = parse(R"(
[system]
spins = 3
[field]
hx = [1.0, 2.0, 3.0]
[measurement]
observables = ["norm"]
[[step]]
kind = "measure"
)");
  REQUIRE(cfg.baseFields[0].x == 1.0);
  REQUIRE(cfg.baseFields[1].x == 2.0);
  REQUIRE(cfg.baseFields[2].x == 3.0);
}

TEST_CASE("a field array of the wrong length is rejected", "[config][errors]") {
  REQUIRE_THROWS_AS(parse(R"(
[system]
spins = 3
[field]
hx = [1.0, 2.0]
[measurement]
observables = ["norm"]
[[step]]
kind = "measure"
)"),
                    ConfigError);
}

TEST_CASE("repeat steps nest", "[config][sequence]") {
  // The legacy format refused this outright: chparsgenPDDG.f returns -987 for
  // a nested @CYCLE, because the cycle was a pair of jump targets rather than
  // a structure.
  const RunConfig cfg = parse(R"(
[system]
spins = 2
[measurement]
observables = ["norm"]

[[step]]
kind = "repeat"
count = 3
body = [
  { kind = "repeat", count = 4, body = [ { kind = "measure" } ] },
]
)");
  REQUIRE(cfg.sequence.measurementCount() == 12);
}

TEST_CASE("a step's field override falls back to the baseline", "[config]") {
  // Naming only hx must not silently zero hz.
  const RunConfig cfg = parse(R"(
[system]
spins = 2
[field]
hz = 5.0
[measurement]
observables = ["norm"]
[[step]]
kind = "evolve"
tau = 0.1
hx = 100.0
[[step]]
kind = "measure"
)");
  const auto& step = std::get<Evolve>(cfg.sequence.steps()[0]);
  REQUIRE(step.fields.size() == 2);
  REQUIRE(step.fields[0].x == 100.0);
  REQUIRE(step.fields[0].z == 5.0);
}

TEST_CASE("a sequence that records nothing is rejected", "[config][errors]") {
  REQUIRE_THROWS_AS(parse(R"(
[system]
spins = 2
[measurement]
observables = ["norm"]
[[step]]
kind = "evolve"
tau = 0.1
)"),
                    ConfigError);
}

TEST_CASE("configuration errors name the offending key", "[config][errors]") {
  try {
    parse(R"(
[system]
spins = 2
[measurement]
observables = ["corr:w"]
[[step]]
kind = "measure"
)");
    FAIL("expected a ConfigError");
  } catch (const ConfigError& e) {
    const std::string message = e.what();
    INFO(message);
    REQUIRE(message.find("observables") != std::string::npos);
    REQUIRE(message.find("<test>") != std::string::npos);
  }
}

TEST_CASE("unknown step kinds are reported with the alternatives",
          "[config][errors]") {
  try {
    parse(R"(
[system]
spins = 2
[measurement]
observables = ["norm"]
[[step]]
kind = "teleport"
)");
    FAIL("expected a ConfigError");
  } catch (const ConfigError& e) {
    const std::string message = e.what();
    REQUIRE(message.find("teleport") != std::string::npos);
    REQUIRE(message.find("evolve") != std::string::npos);
  }
}

TEST_CASE("malformed TOML is reported as a syntax error", "[config][errors]") {
  REQUIRE_THROWS_AS(parse("[system\nspins = 2"), ConfigError);
}

TEST_CASE("too many spins are rejected", "[config][errors]") {
  REQUIRE_THROWS_AS(parse(R"(
[system]
spins = 99
[measurement]
observables = ["norm"]
[[step]]
kind = "measure"
)"),
                    ConfigError);
}

TEST_CASE("a realisation is a pure function of config and index",
          "[config][experiment][determinism]") {
  RunConfig cfg = parse(kMinimal);
  cfg.seed = 4242;

  VectorRecorder first;
  VectorRecorder second;
  runRealization(cfg, 9, first);
  runRealization(cfg, 9, second);

  REQUIRE(first.values() == second.values());
  REQUIRE(first.times() == second.times());

  VectorRecorder other;
  runRealization(cfg, 10, other);
  REQUIRE(other.values() != first.values());
}

TEST_CASE("results carry self-describing column headings",
          "[config][experiment]") {
  RunConfig cfg = parse(kMinimal);
  cfg.realizations = 4;

  LocalRunner runner(2);
  const EnsembleResult result = runExperiment(cfg, runner);

  REQUIRE(result.columns == std::vector<std::string>{"corr_z_re", "corr_z_im"});
  REQUIRE(result.times.size() == 2);
  REQUIRE(result.statistics.count() == 4);

  std::ostringstream csv;
  writeCsv(result, csv);
  const std::string text = csv.str();
  REQUIRE(text.starts_with("time,corr_z_re_mean,corr_z_re_sd,"
                           "corr_z_re_stderr,"));
}

TEST_CASE("the correlator starts at the spin count", "[config][experiment]") {
  // At t = 0 the correlator is <A psi|A psi> = ||A psi||^2, which for
  // A = sum sigma_z on a random state averages to the number of spins.
  RunConfig cfg = parse(R"(
[system]
spins = 6
[measurement]
reference_axis = "z"
observables = ["corr:z"]
[ensemble]
realizations = 32
seed = 5150
dipolar_geometry = false
[[step]]
kind = "measure"
)");

  LocalRunner runner(4);
  const EnsembleResult result = runExperiment(cfg, runner);
  REQUIRE(result.statistics.mean()[0] == Catch::Approx(6.0).margin(0.5));
}
