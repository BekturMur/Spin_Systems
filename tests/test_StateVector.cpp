#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "DenseReference.hpp"
#include "spinsim/core/StateVector.hpp"

using namespace spinsim;
using namespace spinsim::test;

TEST_CASE("a fresh state has 2^L amplitudes, all zero", "[state]") {
  const StateVector psi(5);
  REQUIRE(psi.nspins() == 5);
  REQUIRE(psi.size() == 32);
  for (std::size_t k = 0; k < psi.size(); ++k) {
    REQUIRE(psi.re(k) == 0.0);
    REQUIRE(psi.im(k) == 0.0);
  }
}

TEST_CASE("amplitude storage is interleaved and 64-byte aligned", "[state]") {
  const StateVector psi(8);
  REQUIRE(reinterpret_cast<std::uintptr_t>(psi.data()) % 64 == 0);
  // Amplitude k lives at 2k and 2k+1, so the backing array is twice as long
  // as the state space. This is what lets the flip-flop kernel move a whole
  // amplitude as one 128-bit quantity.
  REQUIRE(psi.amplitudes().size() == 2 * psi.size());
  // Amplitude 3 must occupy data()[6] and data()[7], adjacently.
  StateVector probe(8);
  probe.re(3) = 1.25;
  probe.im(3) = -4.5;
  REQUIRE(probe.data()[6] == 1.25);
  REQUIRE(probe.data()[7] == -4.5);
}

TEST_CASE("setBasisState collapses onto one amplitude", "[state]") {
  StateVector psi(3);
  psi.setBasisState(5);
  REQUIRE(psi.re(5) == 1.0);
  REQUIRE(psi.norm() == Catch::Approx(1.0));
  REQUIRE_THROWS_AS(psi.setBasisState(8), std::out_of_range);
}

TEST_CASE("norm and normalize agree", "[state]") {
  StateVector psi(4);
  psi.re(0) = 3.0;
  psi.im(1) = 4.0;
  REQUIRE(psi.norm() == Catch::Approx(5.0));

  REQUIRE(psi.normalize() == Catch::Approx(5.0));
  REQUIRE(psi.norm() == Catch::Approx(1.0));

  StateVector zero(2);
  REQUIRE_THROWS_AS(zero.normalize(), std::domain_error);
}

TEST_CASE("innerProduct computes <other|this>", "[state]") {
  // <other|this> with other = |0> + i|1>, this = 2|0> - i|1>:
  // conj(1)*2 + conj(i)*(-i) = 2 + (-i)(-i) = 2 - 1 = 1.
  StateVector self(1);
  self.re(0) = 2.0;
  self.im(1) = -1.0;

  StateVector other(1);
  other.re(0) = 1.0;
  other.im(1) = 1.0;

  const auto [re, im] = self.innerProduct(other);
  REQUIRE(re == Catch::Approx(1.0));
  REQUIRE(im == Catch::Approx(0.0).margin(1e-15));
}

TEST_CASE("innerProduct with itself is the squared norm", "[state]") {
  const StateVector psi = randomState(5, 0x2024);
  const auto [re, im] = psi.innerProduct(psi);
  REQUIRE(re == Catch::Approx(1.0));
  REQUIRE(im == Catch::Approx(0.0).margin(1e-14));
}

TEST_CASE("memory scales with the run, not with kMaxSpins", "[state]") {
  // The legacy code sized every local buffer at 2^24 regardless of the run;
  // a ten-spin state here is 1024 amplitudes, not 16.7 million.
  const StateVector small(10);
  REQUIRE(small.size() == 1024);
}
