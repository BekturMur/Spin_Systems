#include <charconv>
#include <chrono>
#include <exception>
#include <filesystem>
#include <format>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "spinsim/io/Config.hpp"
#include "spinsim/io/Experiment.hpp"

namespace {

using namespace spinsim;

constexpr std::string_view kUsage = R"(spinsim - Chebyshev spin dynamics

Usage:
  spinsim run <config.toml> [--threads N] [--out FILE] [--realizations N]
  spinsim validate <config.toml>
  spinsim describe <config.toml>

Options:
  --threads N       worker threads; 0 or omitted uses every core
  --out FILE        write CSV here instead of standard output
  --realizations N  override ensemble.realizations from the config

Results are averaged over the ensemble, with the standard deviation and the
uncertainty of the mean reported alongside each value.
)";

struct Options {
  std::string command;
  std::filesystem::path config;
  std::size_t threads = 0;
  std::optional<std::size_t> realizations;
  std::optional<std::filesystem::path> output;
};

std::size_t parseCount(std::string_view text, std::string_view flag) {
  std::size_t value = 0;
  const auto* last = text.data() + text.size();
  const auto [ptr, ec] = std::from_chars(text.data(), last, value);
  if (ec != std::errc{} || ptr != last) {
    throw std::invalid_argument(
        std::format("{} expects a non-negative integer, got '{}'", flag, text));
  }
  return value;
}

Options parseArguments(std::span<const std::string_view> args) {
  if (args.size() < 2) throw std::invalid_argument("no command given");

  Options options;
  options.command = args[0];
  options.config = args[1];

  for (std::size_t i = 2; i < args.size(); ++i) {
    const std::string_view flag = args[i];
    const bool needsValue =
        flag == "--threads" || flag == "--out" || flag == "--realizations";
    if (needsValue && i + 1 >= args.size()) {
      throw std::invalid_argument(std::format("{} needs a value", flag));
    }
    if (flag == "--threads") {
      options.threads = parseCount(args[++i], flag);
    } else if (flag == "--realizations") {
      options.realizations = parseCount(args[++i], flag);
    } else if (flag == "--out") {
      options.output = std::filesystem::path(args[++i]);
    } else {
      throw std::invalid_argument(std::format("unknown option '{}'", flag));
    }
  }
  return options;
}

void describe(const RunConfig& config) {
  std::cout << std::format("spins           {}\n", config.nspins);
  std::cout << std::format("state space     {} amplitudes\n",
                           stateCount(config.nspins));
  std::cout << std::format(
      "initial state   {}\n",
      config.initialState == InitialStateKind::Random ? "random" :
      config.initialState == InitialStateKind::AllUp ? "up" :
      config.initialState == InitialStateKind::AllDown ? "down"
                                                       : "basis");
  std::cout << std::format("realizations    {}\n", config.realizations);
  std::cout << std::format("seed            {}\n", config.seed);
  std::cout << std::format(
      "geometry        {}\n",
      config.generateGeometry ? "dipolar 2D, per realization" : "from config");
  std::cout << std::format("field disorder  {:g}\n", config.fieldDisorder);
  std::cout << std::format("reference axis  {}\n",
                           axisName(config.referenceAxis));
  std::cout << std::format("epsilon         {:g}\n", config.epsilon);
  std::cout << std::format("measurements    {} per realization\n",
                           config.sequence.measurementCount());
  std::cout << "observables    ";
  for (const std::string& spec : config.observables) std::cout << ' ' << spec;
  std::cout << '\n';
}

int run(const Options& options) {
  RunConfig config = loadConfig(options.config);
  if (options.realizations) {
    config.realizations = *options.realizations;
    config.validate(options.config.string());
  }

  if (options.command == "validate") {
    std::cout << std::format("{}: ok\n", options.config.string());
    return 0;
  }
  if (options.command == "describe") {
    describe(config);
    return 0;
  }
  if (options.command != "run") {
    throw std::invalid_argument(
        std::format("unknown command '{}'", options.command));
  }

  LocalRunner runner(options.threads);
  const auto started = std::chrono::steady_clock::now();
  const EnsembleResult result = runExperiment(config, runner);
  const auto elapsed = std::chrono::steady_clock::now() - started;

  if (options.output) {
    writeCsv(result, *options.output);
    std::cerr << std::format(
        "{} realizations on {} threads in {:.2f}s -> {}\n",
        config.realizations, runner.threads(),
        std::chrono::duration<double>(elapsed).count(), options.output->string());
  } else {
    writeCsv(result, std::cout);
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  const std::vector<std::string_view> args(argv + 1, argv + argc);
  if (args.empty() || args[0] == "-h" || args[0] == "--help") {
    std::cout << kUsage;
    return args.empty() ? 2 : 0;
  }

  try {
    return run(parseArguments(args));
  } catch (const std::exception& e) {
    std::cerr << std::format("spinsim: {}\n", e.what());
    return 1;
  }
}
