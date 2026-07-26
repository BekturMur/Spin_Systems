#include "spinsim/io/Config.hpp"

#include <format>

#include <toml++/toml.hpp>

namespace spinsim {
namespace {

/// Reads a required value, reporting the key by name when it is missing or of
/// the wrong type.
template <typename T, typename View>
T require(const View& node, std::string_view key, std::string_view origin) {
  const std::optional<T> value = node[key].template value<T>();
  if (!value) {
    throw ConfigError(origin, key, "missing, or not of the expected type");
  }
  return *value;
}

template <typename T, typename View>
T optional(const View& node, std::string_view key, T fallback) {
  std::optional<T> value = node[key].template value<T>();
  return value ? *std::move(value) : std::move(fallback);
}

/// A field component may be given as one number for every spin, or as an array
/// with one entry per spin.
template <typename View>
void readComponent(const View& node, std::string_view key, std::size_t nspins,
                   std::string_view origin, std::vector<double>& out) {
  out.assign(nspins, 0.0);
  const auto& entry = node[key];
  if (!entry) return;

  if (const std::optional<double> scalar = entry.template value<double>()) {
    out.assign(nspins, *scalar);
    return;
  }
  const toml::array* list = entry.as_array();
  if (list == nullptr) {
    throw ConfigError(origin, key, "expected a number or an array of numbers");
  }
  if (list->size() != nspins) {
    throw ConfigError(origin, key,
                      std::format("has {} entries but there are {} spins",
                                  list->size(), nspins));
  }
  for (std::size_t i = 0; i < nspins; ++i) {
    const std::optional<double> v = (*list)[i].value<double>();
    if (!v) throw ConfigError(origin, key, "contains a non-numeric entry");
    out[i] = *v;
  }
}

template <typename View>
std::vector<Field> readFields(const View& node, std::size_t nspins,
                              std::string_view origin) {
  std::vector<double> hx;
  std::vector<double> hy;
  std::vector<double> hz;
  readComponent(node, "hx", nspins, origin, hx);
  readComponent(node, "hy", nspins, origin, hy);
  readComponent(node, "hz", nspins, origin, hz);

  std::vector<Field> fields(nspins);
  for (std::size_t s = 0; s < nspins; ++s) {
    fields[s] = Field{hx[s], hy[s], hz[s]};
  }
  return fields;
}

/// True when the table names any field component, i.e. wants to override the
/// system's baseline for the duration of a step.
bool mentionsField(const toml::table& step) {
  return step.contains("hx") || step.contains("hy") || step.contains("hz");
}

std::vector<Step> readSteps(const toml::array& list, std::size_t nspins,
                            const std::vector<Field>& baseFields,
                            std::string_view origin, std::size_t depth);

Step readStep(const toml::table& table, std::size_t nspins,
              const std::vector<Field>& baseFields, std::string_view origin,
              std::size_t depth) {
  const auto view = toml::node_view<const toml::node>(table);
  const std::string kind = require<std::string>(view, "kind", origin);

  if (kind == "measure") return Measure{};

  if (kind == "evolve") {
    Evolve step;
    step.tau = require<double>(view, "tau", origin);
    const auto steps = optional<std::int64_t>(view, "steps", 1);
    if (steps < 0) throw ConfigError(origin, "steps", "must not be negative");
    step.nsteps = static_cast<std::size_t>(steps);

    const std::string mode = optional<std::string>(view, "mode", "real");
    if (mode == "real") {
      step.mode = TimeMode::Real;
    } else if (mode == "imaginary") {
      step.mode = TimeMode::Imaginary;
    } else {
      throw ConfigError(origin, "mode", "expected \"real\" or \"imaginary\"");
    }
    step.autonormalize = optional<bool>(view, "autonormalize", false);

    if (mentionsField(table)) {
      // Components not named keep the system's baseline value.
      std::vector<double> hx;
      std::vector<double> hy;
      std::vector<double> hz;
      readComponent(view, "hx", nspins, origin, hx);
      readComponent(view, "hy", nspins, origin, hy);
      readComponent(view, "hz", nspins, origin, hz);
      step.fields.resize(nspins);
      for (std::size_t s = 0; s < nspins; ++s) {
        step.fields[s] = Field{
            table.contains("hx") ? hx[s] : baseFields[s].x,
            table.contains("hy") ? hy[s] : baseFields[s].y,
            table.contains("hz") ? hz[s] : baseFields[s].z,
        };
      }
    }
    return step;
  }

  if (kind == "pulse") {
    const std::string spec = require<std::string>(view, "rotation", origin);
    try {
      return Pulse{makeRotation(spec, nspins)};
    } catch (const std::invalid_argument& e) {
      throw ConfigError(origin, "rotation", e.what());
    }
  }

  if (kind == "repeat") {
    if (depth > 32) throw ConfigError(origin, "repeat", "nested too deeply");
    Repeat repeat;
    const auto count = require<std::int64_t>(view, "count", origin);
    if (count < 0) throw ConfigError(origin, "count", "must not be negative");
    repeat.count = static_cast<std::size_t>(count);

    const toml::array* body = table["body"].as_array();
    if (body == nullptr) {
      throw ConfigError(origin, "body", "a repeat needs an array of steps");
    }
    repeat.body = readSteps(*body, nspins, baseFields, origin, depth + 1);
    return repeat;
  }

  throw ConfigError(origin, "kind",
                    std::format("unknown step kind '{}'; expected evolve, "
                                "pulse, measure or repeat",
                                kind));
}

std::vector<Step> readSteps(const toml::array& list, std::size_t nspins,
                            const std::vector<Field>& baseFields,
                            std::string_view origin, std::size_t depth) {
  std::vector<Step> steps;
  steps.reserve(list.size());
  for (const toml::node& node : list) {
    const toml::table* table = node.as_table();
    if (table == nullptr) {
      throw ConfigError(origin, "step", "every step must be a table");
    }
    steps.push_back(readStep(*table, nspins, baseFields, origin, depth));
  }
  return steps;
}

}  // namespace

ConfigError::ConfigError(std::string_view file, std::string_view key,
                         std::string_view problem)
    : std::runtime_error(
          std::format("{}: '{}' {}", file, key, problem)) {}

void RunConfig::validate(std::string_view origin) const {
  if (nspins == 0 || nspins > kMaxSpins) {
    throw ConfigError(origin, "system.spins",
                      std::format("must be between 1 and {}", kMaxSpins));
  }
  if (baseFields.size() != nspins) {
    throw ConfigError(origin, "field", "does not cover every spin");
  }
  if (realizations == 0) {
    throw ConfigError(origin, "ensemble.realizations", "must be at least 1");
  }
  if (!(epsilon > 0.0)) {
    throw ConfigError(origin, "numerics.epsilon", "must be positive");
  }
  if (observables.empty()) {
    throw ConfigError(origin, "measurement.observables", "must not be empty");
  }
  if (sequence.measurementCount() == 0) {
    throw ConfigError(origin, "step",
                      "the sequence records nothing; add a measure step");
  }
  if (initialState == InitialStateKind::BasisState &&
      basisIndex >= stateCount(nspins)) {
    throw ConfigError(origin, "system.basis_index",
                      "outside the 2^L state space");
  }
  for (const std::string& spec : observables) {
    try {
      const auto probe = makeObservable(spec);
      (void)probe;
    } catch (const std::invalid_argument& e) {
      throw ConfigError(origin, "measurement.observables", e.what());
    }
  }
}

RunConfig parseConfig(std::string_view text, std::string_view origin) {
  toml::table root;
  try {
    root = toml::parse(text, origin);
  } catch (const toml::parse_error& e) {
    throw ConfigError(origin, "<syntax>", e.description());
  }

  const toml::table& doc = root;

  RunConfig cfg;
  const auto system = doc["system"];
  const auto spins = require<std::int64_t>(system, "spins", origin);
  if (spins <= 0) throw ConfigError(origin, "system.spins", "must be positive");
  cfg.nspins = static_cast<std::size_t>(spins);

  try {
    cfg.initialState = parseInitialStateKind(
        optional<std::string>(system, "initial_state", "random"));
  } catch (const std::invalid_argument& e) {
    throw ConfigError(origin, "system.initial_state", e.what());
  }
  cfg.basisIndex =
      static_cast<std::size_t>(optional<std::int64_t>(system, "basis_index", 0));

  cfg.baseFields = readFields(doc["field"], cfg.nspins, origin);

  const auto ensemble = doc["ensemble"];
  cfg.realizations = static_cast<std::size_t>(
      optional<std::int64_t>(ensemble, "realizations", 1));
  cfg.seed =
      static_cast<std::uint64_t>(optional<std::int64_t>(ensemble, "seed", 0));
  cfg.generateGeometry = optional<bool>(ensemble, "dipolar_geometry", true);
  cfg.dipolarScale = optional<double>(ensemble, "dipolar_scale", 1.0);
  cfg.fieldDisorder = optional<double>(ensemble, "field_disorder", 0.0);

  const auto measurement = doc["measurement"];
  try {
    cfg.referenceAxis =
        parseAxis(optional<std::string>(measurement, "reference_axis", "z"));
  } catch (const std::invalid_argument& e) {
    throw ConfigError(origin, "measurement.reference_axis", e.what());
  }
  if (const toml::array* list = measurement["observables"].as_array()) {
    for (const toml::node& node : *list) {
      const std::optional<std::string> spec = node.value<std::string>();
      if (!spec) {
        throw ConfigError(origin, "measurement.observables",
                          "entries must be strings");
      }
      cfg.observables.push_back(*spec);
    }
  }

  cfg.epsilon = optional<double>(doc["numerics"], "epsilon", 1.0e-7);

  const toml::array* steps = doc["step"].as_array();
  if (steps == nullptr) {
    throw ConfigError(origin, "step", "missing; a run needs a [[step]] list");
  }
  cfg.sequence =
      Sequence(readSteps(*steps, cfg.nspins, cfg.baseFields, origin, 0));

  cfg.validate(origin);
  return cfg;
}

RunConfig loadConfig(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw ConfigError(path.string(), "<file>", "cannot be opened for reading");
  }
  const std::string text((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());
  return parseConfig(text, path.string());
}

}  // namespace spinsim
