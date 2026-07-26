#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "DenseReference.hpp"
#include "spinsim/core/Hamiltonian.hpp"

using namespace spinsim;
using namespace spinsim::test;

namespace {

constexpr double kTol = 1e-13;

/// A deliberately messy Hamiltonian: every field component populated, couplings
/// that are neither uniform nor nearest-neighbour only.
Hamiltonian sampleHamiltonian(std::size_t nspins) {
  Hamiltonian h(nspins);
  for (SpinIndex s = 0; s < nspins; ++s) {
    const double k = static_cast<double>(s) + 1.0;
    h.setField(s, Field{0.3 * k, -0.2 * k, 1.1 / k});
  }
  for (SpinIndex a = 0; a < nspins; ++a) {
    for (SpinIndex b = a + 1; b < nspins; ++b) {
      const double r = static_cast<double>(b - a);
      h.setCoupling(SpinPair{a, b}, 0.4 / r, -0.15 / r, 0.9 / (r * r));
    }
  }
  return h;
}

}  // namespace

TEST_CASE("applyPhysical matches the dense Hamiltonian", "[hamiltonian]") {
  for (std::size_t nspins = 1; nspins <= 5; ++nspins) {
    const Hamiltonian h = sampleHamiltonian(nspins);
    const StateVector psi = randomState(nspins, 0x4242 + nspins);

    StateVector got(nspins);
    h.applyPhysical(psi, got);

    INFO("nspins=" << nspins);
    REQUIRE(maxAbsDiff(got, denseApply(denseHamiltonian(h), psi)) < kTol);
  }
}

TEST_CASE("applyPhysical overwrites any prior contents", "[hamiltonian]") {
  const Hamiltonian h = sampleHamiltonian(3);
  const StateVector psi = randomState(3, 0x77);

  StateVector fresh(3);
  StateVector dirty = randomState(3, 0x88);
  h.applyPhysical(psi, fresh);
  h.applyPhysical(psi, dirty);

  REQUIRE(maxAbsDiff(fresh, dirty) == 0.0);
}

TEST_CASE("spectralBound really bounds the spectrum", "[hamiltonian]") {
  // The Chebyshev expansion converges only if the rescaled operator has its
  // spectrum inside [-1, 1]. Checked here via the Rayleigh quotient over many
  // random states, which is a lower bound on the true spectral radius.
  for (std::size_t nspins = 1; nspins <= 5; ++nspins) {
    const Hamiltonian h = sampleHamiltonian(nspins);
    const double bound = h.spectralBound();
    REQUIRE(bound > 0.0);

    StateVector hpsi(nspins);
    for (std::uint64_t trial = 0; trial < 32; ++trial) {
      const StateVector psi = randomState(nspins, 0x900 + trial * 7 + nspins);
      h.applyPhysical(psi, hpsi);
      // |<psi|H|psi>| <= ||H psi|| <= bound for a normalised psi.
      INFO("nspins=" << nspins << " trial=" << trial);
      REQUIRE(hpsi.norm() <= bound * (1.0 + 1e-12));
    }
  }
}

TEST_CASE("applyTwiceNormalized is -2H/bound", "[hamiltonian]") {
  // The sign and the factor of two are inherited from the legacy code
  // (hefx = -Hx/emax with the comment "= -0.5*H*2"); this pins both down so a
  // future refactor cannot quietly flip them.
  for (std::size_t nspins = 1; nspins <= 4; ++nspins) {
    const Hamiltonian h = sampleHamiltonian(nspins);
    const double bound = h.spectralBound();
    const StateVector psi = randomState(nspins, 0xABC + nspins);

    StateVector scaled(nspins);
    StateVector physical(nspins);
    h.applyTwiceNormalized(bound, psi, scaled);
    h.applyPhysical(psi, physical);

    for (std::size_t k = 0; k < psi.size(); ++k) {
      const double factor = -2.0 / bound;
      INFO("nspins=" << nspins << " k=" << k);
      REQUIRE(scaled.re(k) == Catch::Approx(factor * physical.re(k))
                                    .margin(kTol));
      REQUIRE(scaled.im(k) == Catch::Approx(factor * physical.im(k))
                                    .margin(kTol));
    }
  }
}

TEST_CASE("normalized operator has spectral radius at most 2",
          "[hamiltonian]") {
  for (std::size_t nspins = 1; nspins <= 5; ++nspins) {
    const Hamiltonian h = sampleHamiltonian(nspins);
    const double bound = h.spectralBound();
    StateVector out(nspins);
    for (std::uint64_t trial = 0; trial < 16; ++trial) {
      const StateVector psi = randomState(nspins, 0xF00 + trial + nspins);
      h.applyTwiceNormalized(bound, psi, out);
      INFO("nspins=" << nspins << " trial=" << trial);
      REQUIRE(out.norm() <= 2.0 * (1.0 + 1e-12));
    }
  }
}

TEST_CASE("repeating a coupling updates it in place", "[hamiltonian]") {
  // Matches chparsfPDDG.f, where a second `i,j,...` line for the same pair
  // overwrites the first rather than adding a duplicate term.
  Hamiltonian h(3);
  h.setCoupling(SpinPair{0, 1}, 1.0, 1.0, 1.0);
  h.setCoupling(SpinPair{1, 0}, 2.0, 3.0, 4.0);

  REQUIRE(h.couplings().size() == 1);
  REQUIRE(h.couplings()[0].jx == 2.0);
  REQUIRE(h.couplings()[0].jy == 3.0);
  REQUIRE(h.couplings()[0].jz == 4.0);
}

TEST_CASE("out-of-range spins are rejected", "[hamiltonian]") {
  Hamiltonian h(3);
  REQUIRE_THROWS_AS(h.setField(3, Field{}), std::out_of_range);
  REQUIRE_THROWS_AS(h.setCoupling(SpinPair{0, 3}, 1, 1, 1), std::out_of_range);
  REQUIRE_THROWS_AS(Hamiltonian(0), std::out_of_range);
}
