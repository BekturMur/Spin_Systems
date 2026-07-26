#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>
#include <numbers>

#include "DenseReference.hpp"
#include "spinsim/core/Propagator.hpp"
#include "spinsim/core/SpinOperators.hpp"

using namespace spinsim;
using namespace spinsim::test;

namespace {

/// Interacting, anisotropic, and with every field component populated -- a
/// Hamiltonian with no accidental symmetry to hide a sign error behind.
Hamiltonian sampleHamiltonian(std::size_t nspins) {
  Hamiltonian h(nspins);
  for (SpinIndex s = 0; s < nspins; ++s) {
    const double k = static_cast<double>(s) + 1.0;
    h.setField(s, Field{0.4 * k, -0.3 / k, 0.9});
  }
  for (SpinIndex a = 0; a < nspins; ++a) {
    for (SpinIndex b = a + 1; b < nspins; ++b) {
      const double r = static_cast<double>(b - a);
      h.setCoupling(SpinPair{a, b}, 0.5 / r, -0.2 / r, 0.8 / (r * r));
    }
  }
  return h;
}

}  // namespace

TEST_CASE("real-time propagation matches the exact matrix exponential",
          "[propagator][exact]") {
  // The decisive check of the whole port: an independent algorithm (dense
  // eigendecomposition) for the same operator.
  //
  // The sign is the one the legacy code produces. chstepsPDDGnmr.f:175 sets
  // hef = -H/emax, so the expansion sums (-i)^n T_n(-H/emax) = (+i)^n
  // T_n(H/emax), giving exp(+i H t). The legacy convention therefore evolves
  // with the opposite sign to the textbook Schroedinger equation -- equivalently
  // it treats the tabulated fields as -H. Every observable in this code is
  // built from |<a|b>| or from correlators symmetric under time reversal, which
  // is why the difference was never visible in the published results.
  const double epsilon = 1e-12;

  for (std::size_t nspins = 1; nspins <= 6; ++nspins) {
    const Hamiltonian h = sampleHamiltonian(nspins);
    const StateVector psi0 = randomState(nspins, 0x9001 + nspins);
    const DenseMatrix dense = denseHamiltonian(h);

    for (const double tau : {0.05, 0.3, 1.0}) {
      StateVector psi = psi0;
      Propagator prop(nspins, epsilon);
      prop.evolve(h, TimeMode::Real, psi, tau, 1);

      const StateVector want =
          denseExpApply(dense, Complex{0.0, 1.0} * tau, psi0);

      INFO("nspins=" << nspins << " tau=" << tau
                     << " order=" << prop.lastOrder());
      REQUIRE(maxAbsDiff(psi, want) < 1e-10);
    }
  }
}

TEST_CASE("many small steps agree with one exact long step",
          "[propagator][exact]") {
  constexpr std::size_t kSpins = 4;
  const Hamiltonian h = sampleHamiltonian(kSpins);
  const StateVector psi0 = randomState(kSpins, 0x5150);

  StateVector psi = psi0;
  Propagator prop(kSpins, 1e-12);
  prop.evolve(h, TimeMode::Real, psi, 0.01, 200);

  const StateVector want =
      denseExpApply(denseHamiltonian(h), Complex{0.0, 1.0} * 2.0, psi0);
  REQUIRE(maxAbsDiff(psi, want) < 1e-9);
}

TEST_CASE("a negative step undoes a positive one", "[propagator]") {
  // Exercises the J_n(-a) = (-1)^n J_n(a) sign flip on the odd coefficients.
  constexpr std::size_t kSpins = 4;
  const Hamiltonian h = sampleHamiltonian(kSpins);
  const StateVector psi0 = randomState(kSpins, 0x2718);

  StateVector psi = psi0;
  Propagator prop(kSpins, 1e-12);
  prop.evolve(h, TimeMode::Real, psi, 0.7, 1);
  REQUIRE(maxAbsDiff(psi, psi0) > 1e-3);  // it really did move

  prop.evolve(h, TimeMode::Real, psi, -0.7, 1);
  REQUIRE(maxAbsDiff(psi, psi0) < 1e-11);
}

TEST_CASE("real-time propagation conserves norm and energy",
          "[propagator][conservation]") {
  // Without autonormalisation, so this measures the expansion's accuracy rather
  // than hiding it. The legacy `autonorm` switch could mask exactly this.
  constexpr std::size_t kSpins = 5;
  const Hamiltonian h = sampleHamiltonian(kSpins);
  StateVector psi = randomState(kSpins, 0x3141);

  const auto energy = [&](const StateVector& s) {
    StateVector hs(kSpins);
    h.applyPhysical(s, hs);
    return s.innerProduct(hs).first;
  };

  const double norm0 = psi.norm();
  const double energy0 = energy(psi);

  Propagator prop(kSpins, 1e-12);
  for (int block = 0; block < 100; ++block) {
    prop.evolve(h, TimeMode::Real, psi, 0.02, 100);
  }

  REQUIRE(psi.norm() == Catch::Approx(norm0).margin(1e-11));
  REQUIRE(energy(psi) == Catch::Approx(energy0).margin(1e-10));
}

TEST_CASE("a single spin in a transverse field precesses", "[propagator]") {
  // H = Hx Sx with S = sigma/2, so the field drives rotation at angular
  // frequency Hx. Starting from |down>, <sigma_z> returns -cos(Hx t) under the
  // legacy sign convention, which is cos-like either way.
  Hamiltonian h(1);
  const double hx = 2.0;
  h.setField(0, Field{hx, 0.0, 0.0});

  StateVector psi(1);
  psi.setBasisState(0);  // spin down

  Propagator prop(1, 1e-13);
  const double dt = 0.05;
  double t = 0.0;

  for (int step = 0; step < 60; ++step) {
    prop.evolve(h, TimeMode::Real, psi, dt, 1);
    t += dt;

    StateVector zpsi(1);
    applySigmaZ(psi, 0, zpsi);
    const double sz = psi.innerProduct(zpsi).first;

    INFO("t=" << t);
    REQUIRE(sz == Catch::Approx(-std::cos(hx * t)).margin(1e-9));
  }
}

TEST_CASE("a negligibly short step leaves the state untouched",
          "[propagator]") {
  constexpr std::size_t kSpins = 3;
  const Hamiltonian h = sampleHamiltonian(kSpins);
  const StateVector psi0 = randomState(kSpins, 0x606);

  StateVector psi = psi0;
  Propagator prop(kSpins, 1e-7);
  prop.evolve(h, TimeMode::Real, psi, 1e-12, 1);

  REQUIRE(prop.lastOrder() == 0);
  REQUIRE(maxAbsDiff(psi, psi0) == 0.0);
}

TEST_CASE("imaginary-time propagation matches exp(H t)",
          "[propagator][exact][imaginary]") {
  // chimstepsPDDGnmr.f expands in I_n with alternating signs, which together
  // with Hhat = -H/emax gives exp(+H t) -- amplifying the top of the spectrum.
  constexpr std::size_t kSpins = 4;
  const Hamiltonian h = sampleHamiltonian(kSpins);
  const StateVector psi0 = randomState(kSpins, 0x1618);

  for (const double tau : {0.05, 0.4}) {
    StateVector psi = psi0;
    Propagator prop(kSpins, 1e-12);
    prop.evolve(h, TimeMode::Imaginary, psi, tau, 1);

    const StateVector want =
        denseExpApply(denseHamiltonian(h), Complex{tau, 0.0}, psi0);
    INFO("tau=" << tau << " order=" << prop.lastOrder());
    REQUIRE(maxAbsDiff(psi, want) < 1e-10);
  }
}

TEST_CASE("expansion order grows with the step length", "[propagator]") {
  constexpr std::size_t kSpins = 3;
  const Hamiltonian h = sampleHamiltonian(kSpins);
  StateVector psi = randomState(kSpins, 0x4004);

  Propagator prop(kSpins, 1e-10);
  prop.evolve(h, TimeMode::Real, psi, 0.1, 1);
  const std::size_t shortOrder = prop.lastOrder();
  prop.evolve(h, TimeMode::Real, psi, 2.0, 1);
  const std::size_t longOrder = prop.lastOrder();

  REQUIRE(shortOrder > 0);
  REQUIRE(longOrder > shortOrder);
}

TEST_CASE("an over-long step is reported, not silently wrong",
          "[propagator][errors]") {
  // The legacy code returned status -204 and smuggled the diagnostics out
  // through the results array; here the same information rides on the
  // exception.
  // alpha = bound * |tau| far exceeds any expansion of length maxOrder.
  REQUIRE_THROWS_AS(
      ChebyshevExpansion::build(TimeMode::Real, 1.0, 1e6, 1e-7, 100),
      ChebyshevConvergenceError);

  // The message carries what a user needs to act on, rather than a bare code.
  try {
    ChebyshevExpansion::build(TimeMode::Real, 1.0, 1e6, 1e-7, 100);
  } catch (const ChebyshevConvergenceError& e) {
    REQUIRE(std::string(e.what()).find("tau") != std::string::npos);
    REQUIRE(e.spectralBound() == 1.0);
    REQUIRE(e.epsilon() == 1e-7);
  }
}

TEST_CASE("mismatched state and Hamiltonian sizes are rejected",
          "[propagator][errors]") {
  const Hamiltonian h = sampleHamiltonian(3);
  StateVector psi(4);
  Propagator prop(4);
  REQUIRE_THROWS_AS(prop.evolve(h, TimeMode::Real, psi, 0.1, 1),
                    std::invalid_argument);
}
