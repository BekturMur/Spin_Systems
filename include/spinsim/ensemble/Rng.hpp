#pragma once

#include <cstdint>
#include <random>

namespace spinsim {

/// Random source for one Monte-Carlo realisation.
///
/// This replaces the legacy scheme, which had three separate problems:
///
///  * `srand` was never called anywhere, so `rand()` always began from the
///    default seed.
///  * "Seeding" meant discarding values in a loop -- `do i=1,lociRand;
///    xx=rand()` in gendipnmr2D.f -- costing up to ~2e6 calls per realisation.
///  * The seed accumulated across jobs. chebNMR2D-v6.f:468 does
///    `iRandT(ifile) = iRandT(ifile) + it*11489` on a worker-local copy that is
///    never reset, so a realisation's random numbers depended on how many jobs
///    that particular MPI rank had already handled, and therefore on
///    scheduling. Two runs of the same input could not be expected to agree.
///
/// Here each realisation derives its stream from (base seed, realisation
/// index) alone. Results do not depend on thread count, execution order, or
/// how the work was distributed.
class Rng {
 public:
  explicit Rng(std::uint64_t seed) : engine_(mix(seed)) {}

  /// The stream for one realisation of an ensemble.
  [[nodiscard]] static Rng forRealization(std::uint64_t baseSeed,
                                          std::size_t index) {
    // SplitMix64 finalisation, so that adjacent indices give unrelated streams
    // rather than merely offset ones.
    return Rng(mix(baseSeed ^ (mix(index + 1) * 0x9E3779B97F4A7C15ULL)));
  }

  /// Uniform on [0, 1).
  [[nodiscard]] double uniform() {
    // 53 significant bits, the most a double can carry.
    return static_cast<double>(engine_() >> 11) * 0x1.0p-53;
  }

  /// Standard normal, by the Box-Muller transform the legacy code also used
  /// (chinitPDDGnmr.f:67-71).
  [[nodiscard]] double gaussian();

 private:
  [[nodiscard]] static std::uint64_t mix(std::uint64_t z) {
    z += 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
  }

  std::mt19937_64 engine_;
  double spare_ = 0.0;
  bool hasSpare_ = false;
};

}  // namespace spinsim
