#pragma once

#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include "spinsim/ensemble/Ensemble.hpp"
#include "spinsim/sequence/Sequence.hpp"

namespace spinsim {

/// Thrown for any problem in a configuration file, naming the offending key.
///
/// The legacy parser reported trouble as a negative integer written to a log
/// file -- -11, -83, -1701, -4444 -- and carried on reading. Worse, it matched
/// keywords by substring, so `index(curline,'END')` closed a block on any line
/// containing those three letters anywhere.
class ConfigError : public std::runtime_error {
 public:
  ConfigError(std::string_view file, std::string_view key,
              std::string_view problem);
};

/// A complete description of one experiment.
///
/// This is the object that replaces forty-seven directories of near-identical
/// Fortran. Everything that used to require editing and recompiling the source
/// -- which correlator to measure, how many realisations to average, how many
/// spins -- is a field here.
struct RunConfig {
  // -- system ------------------------------------------------------------
  std::size_t nspins = 0;
  std::vector<Field> baseFields;  ///< one per spin
  InitialStateKind initialState = InitialStateKind::Random;
  std::size_t basisIndex = 0;

  // -- ensemble ----------------------------------------------------------
  std::size_t realizations = 1;
  std::uint64_t seed = 0;
  bool generateGeometry = true;  ///< run the dipolar 2D placement per realisation
  double dipolarScale = 1.0;
  double fieldDisorder = 0.0;  ///< sigma of the Gaussian offset on Hz

  // -- measurement -------------------------------------------------------
  Axis referenceAxis = Axis::Z;
  std::vector<std::string> observables;

  // -- numerics ----------------------------------------------------------
  double epsilon = 1.0e-7;

  // -- program -----------------------------------------------------------
  Sequence sequence{{}};

  /// Checks internal consistency, throwing ConfigError on the first problem.
  void validate(std::string_view origin = "<config>") const;
};

/// Reads a TOML configuration file.
[[nodiscard]] RunConfig loadConfig(const std::filesystem::path& path);

/// Parses TOML held in memory. `origin` only appears in error messages.
[[nodiscard]] RunConfig parseConfig(std::string_view text,
                                    std::string_view origin);

}  // namespace spinsim
