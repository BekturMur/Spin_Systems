// Micro-benchmark for the propagator hot path.
//
// Not a test: it reports timings so optimisation work has a baseline to argue
// against. Built as a separate binary, run by hand.
//
//   ./build/bench_propagator [nspins] [steps]

#include <charconv>
#include <chrono>
#include <cstdio>
#include <string_view>
#include <vector>

#include "spinsim/core/Propagator.hpp"
#include "spinsim/core/SpinOperators.hpp"
#include "spinsim/ensemble/Ensemble.hpp"

using namespace spinsim;

namespace {

double seconds(auto duration) {
  return std::chrono::duration<double>(duration).count();
}

/// Times `body`, reporting the best of a few passes to suppress scheduler noise.
template <typename Body>
double timeBest(int repeats, Body body) {
  double best = 1e300;
  for (int r = 0; r < repeats; ++r) {
    const auto start = std::chrono::steady_clock::now();
    body();
    best = std::min(best, seconds(std::chrono::steady_clock::now() - start));
  }
  return best;
}

std::size_t parse(std::string_view text, std::size_t fallback) {
  std::size_t value = 0;
  const auto* last = text.data() + text.size();
  const auto [ptr, ec] = std::from_chars(text.data(), last, value);
  return (ec == std::errc{} && ptr == last) ? value : fallback;
}

}  // namespace

int main(int argc, char** argv) {
  const std::size_t nspins = argc > 1 ? parse(argv[1], 14) : 14;
  const std::size_t steps = argc > 2 ? parse(argv[2], 20) : 20;

  // The same shape the CPMG configuration produces: all-to-all dipolar
  // couplings from the 2D generator, plus Gaussian disorder on Hz.
  Rng rng = Rng::forRealization(20240726, 0);
  Hamiltonian h(nspins);
  for (const Coupling& c : generateDipolarLattice2D(nspins, rng, 1.0)) {
    h.setCoupling(c.pair, c.jx, c.jy, c.jz);
  }
  applyFieldDisorder(h, rng, 100.0);

  StateVector psi(nspins);
  prepareInitialState(psi, InitialStateKind::Random, rng);

  const double bound = h.spectralBound();
  const std::size_t nstates = psi.size();

  std::printf("spins            %zu  (%zu amplitudes, %.1f MB per state)\n",
              nspins, nstates,
              2.0 * static_cast<double>(nstates) * 8.0 / (1024 * 1024));
  std::printf("couplings        %zu\n", h.couplings().size());
  std::printf("spectral bound   %.4g\n", bound);

  // -- individual kernels ------------------------------------------------
  StateVector in = psi;
  StateVector out(nspins);

  const double tField = timeBest(5, [&] {
    for (SpinIndex s = 0; s < nspins; ++s) {
      axpySigmaX(0.5, in, s, out);
      axpySigmaY(0.5, in, s, out);
      axpySigmaZ(0.5, in, s, out);
    }
  });
  const double tCoupling = timeBest(5, [&] {
    for (const Coupling& c : h.couplings()) {
      axpyCoupling(c.jx, c.jy, c.jz, in, c.pair, out);
    }
  });
  const double tApply = timeBest(5, [&] {
    h.accumulateTwiceNormalized(bound, in, out);
  });

  std::printf("\nper Hamiltonian application:\n");
  std::printf("  all field kernels (3 per spin)  %8.3f ms\n", tField * 1e3);
  std::printf("  all coupling kernels            %8.3f ms\n", tCoupling * 1e3);
  std::printf("  accumulateTwiceNormalized       %8.3f ms\n", tApply * 1e3);
  std::printf("  -> coupling share               %8.1f %%\n",
              100.0 * tCoupling / (tField + tCoupling));

  // -- whole steps -------------------------------------------------------
  const double tau = 0.05;
  Propagator prop(nspins, 1e-7);
  StateVector work = psi;
  prop.evolve(h, TimeMode::Real, work, tau, 1);
  const std::size_t order = prop.lastOrder();

  const double tSteps = timeBest(3, [&] {
    StateVector local = psi;
    prop.evolve(h, TimeMode::Real, local, tau, steps);
  });

  std::printf("\npropagation (tau = %.3g, expansion order %zu):\n", tau, order);
  std::printf("  %zu steps                       %8.1f ms\n", steps,
              tSteps * 1e3);
  std::printf("  per step                        %8.3f ms\n",
              tSteps * 1e3 / static_cast<double>(steps));
  std::printf("  per expansion order             %8.3f ms\n",
              tSteps * 1e3 / static_cast<double>(steps * order));
  std::printf("  amplitude-updates/s             %8.3g\n",
              static_cast<double>(steps * order * nstates) / tSteps);
  return 0;
}
