#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <numbers>

#include "DenseReference.hpp"
#include "spinsim/core/Gate.hpp"
#include "spinsim/core/SpinOperators.hpp"

using namespace spinsim;
using namespace spinsim::test;

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kTol = 1e-13;

}  // namespace

TEST_CASE("rotations preserve the norm", "[gate]") {
  // The check that catches the legacy pi/2 bug. chspecstepsPDDGnmr.f:277 sets
  // ar = sqrt(2) where exp(-i pi/4 sigma_x) needs 1/sqrt(2), so each
  // application there doubles the norm instead of preserving it.
  constexpr std::size_t kSpins = 4;

  for (const double angle : {0.25 * kPi, 0.5 * kPi, kPi, 1.7}) {
    for (const Axis3 axis : {Axis3{1, 0, 0}, Axis3{0, 1, 0}, Axis3{0, 0, 1},
                             Axis3{0.3, -0.5, 0.8}}) {
      StateVector psi = randomState(kSpins, 0xB0A7);
      SpinRotation::onAllSpins(axis, angle, kSpins).apply(psi);
      INFO("angle=" << angle);
      REQUIRE(psi.norm() == Catch::Approx(1.0).margin(1e-14));
    }
  }
}

TEST_CASE("a pi rotation about x is -i sigma_x", "[gate]") {
  // chspecstepsPDDGnmr.f:260 describes exactly this: "pi-X rotation of all
  // spins, incl. factor -I".
  constexpr std::size_t kSpins = 3;
  const StateVector psi0 = randomState(kSpins, 0x9A9A);

  StateVector got = psi0;
  SpinRotation(Axis3{1, 0, 0}, kPi, {1}).apply(got);

  StateVector sx(kSpins);
  applySigmaX(psi0, 1, sx);
  StateVector want(kSpins);
  for (std::size_t k = 0; k < want.size(); ++k) {
    // -i * sx
    want.re(k) = sx.im(k);
    want.im(k) = -sx.re(k);
  }

  REQUIRE(maxAbsDiff(got, want) < kTol);
}

TEST_CASE("two pi rotations return to the start up to a sign", "[gate]") {
  constexpr std::size_t kSpins = 3;
  const StateVector psi0 = randomState(kSpins, 0x1D1D);

  StateVector psi = psi0;
  const SpinRotation gate = SpinRotation::onAllSpins(Axis3{0, 1, 0}, kPi, kSpins);
  gate.apply(psi);
  gate.apply(psi);

  // exp(-i pi n.sigma) = -I on each spin, so three spins give (-1)^3.
  StateVector want = psi0;
  want.negate();
  REQUIRE(maxAbsDiff(psi, want) < kTol);
}

TEST_CASE("a rotation and its inverse cancel", "[gate]") {
  constexpr std::size_t kSpins = 4;
  const StateVector psi0 = randomState(kSpins, 0x3E3E);
  const Axis3 axis{0.6, -0.3, 0.74};

  StateVector psi = psi0;
  SpinRotation::onAllSpins(axis, 0.83, kSpins).apply(psi);
  REQUIRE(maxAbsDiff(psi, psi0) > 1e-3);
  SpinRotation::onAllSpins(axis, -0.83, kSpins).apply(psi);
  REQUIRE(maxAbsDiff(psi, psi0) < kTol);
}

TEST_CASE("a pi-x pulse flips the magnetisation along z", "[gate]") {
  Hamiltonian h(1);
  StateVector psi(1);
  psi.setBasisState(1);  // spin up

  StateVector zpsi(1);
  applySigmaZ(psi, 0, zpsi);
  REQUIRE(psi.innerProduct(zpsi).first == Catch::Approx(1.0));

  SpinRotation(Axis3{1, 0, 0}, kPi, {0}).apply(psi);
  applySigmaZ(psi, 0, zpsi);
  REQUIRE(psi.innerProduct(zpsi).first == Catch::Approx(-1.0));
}

TEST_CASE("rotations act only on their targets", "[gate]") {
  constexpr std::size_t kSpins = 3;
  StateVector psi(kSpins);
  psi.setBasisState(0);  // all down

  // Flip spin 2 (index 1) only; the result must be the basis state with just
  // that bit set, up to the -i phase.
  SpinRotation(Axis3{1, 0, 0}, kPi, {1}).apply(psi);

  REQUIRE(std::abs(psi.im(spinMask(1))) == Catch::Approx(1.0));
  for (std::size_t k = 0; k < psi.size(); ++k) {
    if (k == spinMask(1)) continue;
    REQUIRE(std::hypot(psi.re(k), psi.im(k)) < kTol);
  }
}

TEST_CASE("two-spin gates are checked for unitarity", "[gate]") {
  // The legacy EPR gate at ispecflag 14 swaps the two singlet-triplet
  // components; expressed as a matrix it is a permutation and so unitary.
  TwoSpinGate::Matrix swap{};
  swap[0 * 4 + 0] = 1.0;
  swap[1 * 4 + 2] = 1.0;
  swap[2 * 4 + 1] = 1.0;
  swap[3 * 4 + 3] = 1.0;
  REQUIRE(TwoSpinGate(swap, SpinPair{0, 1}).unitarityDefect() < kTol);

  // The same matrix scaled by sqrt(2) -- the shape of the legacy pi/2 bug --
  // is caught.
  TwoSpinGate::Matrix scaled = swap;
  for (auto& v : scaled) v *= std::numbers::sqrt2;
  REQUIRE(TwoSpinGate(scaled, SpinPair{0, 1}).unitarityDefect() > 0.5);
}

TEST_CASE("a two-spin swap gate permutes amplitudes", "[gate]") {
  TwoSpinGate::Matrix swap{};
  swap[0 * 4 + 0] = 1.0;
  swap[1 * 4 + 2] = 1.0;
  swap[2 * 4 + 1] = 1.0;
  swap[3 * 4 + 3] = 1.0;

  constexpr std::size_t kSpins = 2;
  StateVector psi(kSpins);
  psi.setBasisState(1);  // low spin up

  TwoSpinGate(swap, SpinPair{0, 1}).apply(psi);
  REQUIRE(psi.re(2) == Catch::Approx(1.0));  // now high spin up
  REQUIRE(psi.re(1) == Catch::Approx(0.0).margin(1e-15));
}

TEST_CASE("pulses are built from configuration strings", "[gate]") {
  const SpinRotation all = makeRotation("pi:x", 4);
  REQUIRE(all.targets().size() == 4);
  REQUIRE(all.angle() == Catch::Approx(kPi));

  const SpinRotation one = makeRotation("pi/2:y:3", 4);
  REQUIRE(one.targets() == std::vector<SpinIndex>{2});  // 1-based in config
  REQUIRE(one.angle() == Catch::Approx(0.5 * kPi));
}

TEST_CASE("bad pulse specifications are rejected", "[gate][errors]") {
  REQUIRE_THROWS_AS(makeRotation("pi", 4), std::invalid_argument);
  REQUIRE_THROWS_AS(makeRotation("tau:x", 4), std::invalid_argument);
  REQUIRE_THROWS_AS(makeRotation("pi:w", 4), std::invalid_argument);
  REQUIRE_THROWS_AS(makeRotation("pi:x:0", 4), std::invalid_argument);
  REQUIRE_THROWS_AS(makeRotation("pi:x:9", 4), std::invalid_argument);
  REQUIRE_THROWS_AS(SpinRotation(Axis3{0, 0, 0}, kPi, {0}),
                    std::invalid_argument);
}
