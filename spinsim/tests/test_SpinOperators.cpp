#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "DenseReference.hpp"
#include "spinsim/core/SpinOperators.hpp"

using namespace spinsim;
using namespace spinsim::test;

namespace {

constexpr double kTol = 1e-13;

/// Every (nspins, spin) combination worth checking for L up to 5. Small enough
/// that the dense reference stays cheap, large enough that a spin sits at the
/// bottom, middle and top of the index word.
struct SingleSpinCase {
  std::size_t nspins;
  SpinIndex spin;
};

std::vector<SingleSpinCase> singleSpinCases() {
  std::vector<SingleSpinCase> cases;
  for (std::size_t n = 1; n <= 5; ++n) {
    for (SpinIndex s = 0; s < n; ++s) cases.push_back({n, s});
  }
  return cases;
}

}  // namespace

TEST_CASE("single-spin axpy kernels match dense Pauli matrices",
          "[kernels][operators]") {
  for (const auto& [nspins, spin] : singleSpinCases()) {
    const StateVector psi = randomState(nspins, 0x5eed'0001 + nspins * 31 + spin);
    for (const char axis : {'x', 'y', 'z'}) {
      const double c = 0.75;  // non-unit, to catch a dropped coefficient

      StateVector got(nspins);
      switch (axis) {
        case 'x': axpySigmaX(c, psi, spin, got); break;
        case 'y': axpySigmaY(c, psi, spin, got); break;
        default:  axpySigmaZ(c, psi, spin, got); break;
      }

      const StateVector want =
          denseApply(pauliOperator(nspins, spin, axis) * Complex{c, 0.0}, psi);
      INFO("nspins=" << nspins << " spin=" << spin << " axis=" << axis);
      REQUIRE(maxAbsDiff(got, want) < kTol);
    }
  }
}

TEST_CASE("axpy kernels accumulate rather than overwrite",
          "[kernels][operators]") {
  // The Chebyshev recurrence depends on this: it sums every Hamiltonian term
  // into one output buffer in a single pass.
  constexpr std::size_t kSpins = 4;
  const StateVector psi = randomState(kSpins, 0xA11'CE);

  StateVector fused(kSpins);
  axpySigmaX(1.0, psi, 0, fused);
  axpySigmaZ(1.0, psi, 2, fused);

  DenseMatrix sum = pauliOperator(kSpins, 0, 'x');
  sum += pauliOperator(kSpins, 2, 'z');
  REQUIRE(maxAbsDiff(fused, denseApply(sum, psi)) < kTol);
}

TEST_CASE("apply kernels overwrite the output", "[kernels][operators]") {
  constexpr std::size_t kSpins = 3;
  const StateVector psi = randomState(kSpins, 0xBEEF);

  for (const char axis : {'x', 'y', 'z'}) {
    for (SpinIndex spin = 0; spin < kSpins; ++spin) {
      // Pre-fill with garbage: an overwriting kernel must ignore it.
      StateVector got = randomState(kSpins, 0xDEAD);
      switch (axis) {
        case 'x': applySigmaX(psi, spin, got); break;
        case 'y': applySigmaY(psi, spin, got); break;
        default:  applySigmaZ(psi, spin, got); break;
      }
      INFO("spin=" << spin << " axis=" << axis);
      REQUIRE(maxAbsDiff(got, denseApply(pauliOperator(kSpins, spin, axis),
                                         psi)) < kTol);
    }
  }
}

TEST_CASE("coupling kernel matches dense sigma-sigma products",
          "[kernels][operators][coupling]") {
  // Anisotropic and asymmetric on purpose: equal jx/jy would hide a mix-up
  // between the flip-both and flip-one channels.
  const double jx = 0.3;
  const double jy = -1.1;
  const double jz = 0.7;

  for (std::size_t nspins = 2; nspins <= 5; ++nspins) {
    const StateVector psi = randomState(nspins, 0xC0FFEE + nspins);
    for (SpinIndex a = 0; a < nspins; ++a) {
      for (SpinIndex b = a + 1; b < nspins; ++b) {
        const SpinPair pair{a, b};
        StateVector got(nspins);
        axpyCoupling(jx, jy, jz, psi, pair, got);

        const StateVector want = denseApply(
            couplingOperator(nspins, pair, jx, jy, jz), psi);
        INFO("nspins=" << nspins << " pair=(" << a << "," << b << ")");
        REQUIRE(maxAbsDiff(got, want) < kTol);
      }
    }
  }
}

TEST_CASE("coupling kernel ignores the order spins are named in",
          "[kernels][operators][coupling]") {
  // SpinPair normalises its arguments, so the invariant the legacy code only
  // stated in a comment ("c i0 > j0 !") cannot be violated by a caller.
  constexpr std::size_t kSpins = 4;
  const StateVector psi = randomState(kSpins, 0x1234);

  StateVector forward(kSpins);
  StateVector reversed(kSpins);
  axpyCoupling(0.5, 0.25, -0.75, psi, SpinPair{1, 3}, forward);
  axpyCoupling(0.5, 0.25, -0.75, psi, SpinPair{3, 1}, reversed);

  REQUIRE(maxAbsDiff(forward, reversed) == 0.0);
}

TEST_CASE("kernels satisfy sigma_x sigma_y = i sigma_z",
          "[kernels][operators][algebra]") {
  // Fixes the relative sign of sigma_y without reference to the dense matrices,
  // so a matching error in both would still be caught. Applying x after y means
  // the product sigma_x sigma_y; the result must equal i sigma_z psi.
  constexpr std::size_t kSpins = 3;
  const StateVector psi = randomState(kSpins, 0x51'6D'A1);

  for (SpinIndex spin = 0; spin < kSpins; ++spin) {
    StateVector afterY(kSpins);
    StateVector product(kSpins);
    applySigmaY(psi, spin, afterY);
    applySigmaX(afterY, spin, product);

    StateVector zPsi(kSpins);
    applySigmaZ(psi, spin, zPsi);

    // i * zPsi, written out in split real/imaginary form.
    StateVector want(kSpins);
    for (std::size_t k = 0; k < want.size(); ++k) {
      want.re(k) = -zPsi.im(k);
      want.im(k) = zPsi.re(k);
    }

    INFO("spin=" << spin);
    REQUIRE(maxAbsDiff(product, want) < kTol);
  }
}

TEST_CASE("Pauli kernels are involutions", "[kernels][operators]") {
  // sigma^2 = I for each axis; a two-step round trip must return the input.
  constexpr std::size_t kSpins = 4;
  const StateVector psi = randomState(kSpins, 0x9876);

  for (const char axis : {'x', 'y', 'z'}) {
    for (SpinIndex spin = 0; spin < kSpins; ++spin) {
      StateVector once(kSpins);
      StateVector twice(kSpins);
      switch (axis) {
        case 'x':
          applySigmaX(psi, spin, once);
          applySigmaX(once, spin, twice);
          break;
        case 'y':
          applySigmaY(psi, spin, once);
          applySigmaY(once, spin, twice);
          break;
        default:
          applySigmaZ(psi, spin, once);
          applySigmaZ(once, spin, twice);
          break;
      }
      INFO("spin=" << spin << " axis=" << axis);
      REQUIRE(maxAbsDiff(psi, twice) < kTol);
    }
  }
}
