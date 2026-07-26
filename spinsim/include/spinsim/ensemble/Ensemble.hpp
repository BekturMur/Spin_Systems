#pragma once

#include <cstddef>
#include <span>
#include <string_view>
#include <vector>

#include "spinsim/core/Hamiltonian.hpp"
#include "spinsim/core/Observable.hpp"
#include "spinsim/core/StateVector.hpp"
#include "spinsim/ensemble/Rng.hpp"

namespace spinsim {

/// How the state at t = 0 is prepared.
///
/// Mirrors the legacy `init_state` codes from chparsfPDDG.f. Note that
/// chebsdPDDGnmr.f:96 overrides whatever the .dat file asked for with
/// `init_state = -3` (Random) for the 2D-NMR runs, so Random is the only case
/// those results ever used.
enum class InitialStateKind {
  Random,      ///< legacy -3: a typical state, sampled from the unit sphere
  AllUp,       ///< legacy -1
  AllDown,     ///< legacy -2
  BasisState,  ///< legacy >= 0: a single computational basis state
};

[[nodiscard]] InitialStateKind parseInitialStateKind(std::string_view name);

/// Fills `psi` according to `kind`, consuming randomness only when needed.
void prepareInitialState(StateVector& psi, InitialStateKind kind, Rng& rng,
                         std::size_t basisIndex = 0);

/// Builds the companion state A|psi>, with A = sum_i sigma_axis(i).
///
/// Propagated alongside the main state so observables can form the two-time
/// correlator. chebsdPDDGnmr.f:108-111 does the same, with the axis fixed in
/// source -- `addHSx` in CPMGZ, `addHSz` in Rabi_10_SzSz.
void prepareReferenceState(const StateVector& psi, Axis axis,
                           StateVector& reference);

/// Randomly placed spins in a square, coupled by the dipolar 1/r^3 law.
///
/// Port of gendipnmr2D.f. Spin 1 sits at the origin; the rest are drawn
/// uniformly from a box of half-width 0.5*sqrt(L), rejected if they land within
/// `minSeparation` of an earlier spin in both coordinates at once. That
/// rejection test is the legacy one, and it is looser than a true minimum
/// distance -- two spins may be arbitrarily close diagonally.
///
/// `scale` multiplies the bare 1/r^3 coupling. The transverse components are
/// then fixed at Jx = Jy = -Jz/2, the secular truncation of the dipolar
/// interaction (applied in chebNMR2D-v6.f:483).
[[nodiscard]] std::vector<Coupling> generateDipolarLattice2D(
    std::size_t nspins, Rng& rng, double scale = 1.0,
    double minSeparation = 5.0e-2);

/// Adds a Gaussian random offset to each spin's z field.
///
/// The static disorder that dephases the ensemble; chebNMR2D-v6.f:474-478 draws
/// it with the same Box-Muller construction.
void applyFieldDisorder(Hamiltonian& hamiltonian, Rng& rng, double sigma);

/// Running mean and variance over realisations, one entry per output column.
///
/// The update is Welford's, as in chebNMR2D-v6.f:314-316. What is new is
/// `merge`, which combines two partial accumulators by Chan's formula so that
/// parallel workers can be folded together.
class Accumulator {
 public:
  explicit Accumulator(std::size_t width);

  void add(std::span<const double> sample);
  void merge(const Accumulator& other);

  [[nodiscard]] std::size_t count() const noexcept { return count_; }
  [[nodiscard]] std::size_t width() const noexcept { return mean_.size(); }
  [[nodiscard]] const std::vector<double>& mean() const noexcept {
    return mean_;
  }

  /// Sample standard deviation, sqrt(M2 / (n - 1)). Zero for n < 2.
  [[nodiscard]] std::vector<double> standardDeviation() const;

  /// Uncertainty of the mean, sd / sqrt(n).
  ///
  /// The legacy code wrote sqrt(M2/n) under this heading
  /// (chebNMR2D-v6.f:333, :425), which is neither quantity: it is the
  /// population standard deviation, too small by sqrt(n) to be an error bar on
  /// the average. Both are reported separately here.
  [[nodiscard]] std::vector<double> standardError() const;

 private:
  std::size_t count_ = 0;
  std::vector<double> mean_;
  std::vector<double> m2_;
};

}  // namespace spinsim
