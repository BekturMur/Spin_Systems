#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "DenseReference.hpp"
#include "spinsim/core/Observable.hpp"
#include "spinsim/core/Propagator.hpp"

using namespace spinsim;
using namespace spinsim::test;

namespace {

MeasureContext contextFor(const StateVector& state, const StateVector& reference,
                          const Hamiltonian& h) {
  MeasureContext ctx;
  ctx.state = &state;
  ctx.reference = &reference;
  ctx.hamiltonian = &h;
  return ctx;
}

std::vector<double> record(const Observable& obs, const MeasureContext& ctx) {
  std::vector<double> out(obs.width());
  obs.measure(ctx, out);
  return out;
}

}  // namespace

TEST_CASE("the correlator reproduces the legacy echo expression",
          "[observable]") {
  // chebsdPDDGnmr.f builds phi = sum_i sigma_a psi and takes <pz|phi> via four
  // ddot calls; this checks the same number comes out of the dense reference.
  constexpr std::size_t kSpins = 4;
  const StateVector psi = randomState(kSpins, 0x1111);
  const StateVector ref = randomState(kSpins, 0x2222);
  const Hamiltonian h(kSpins);

  for (const Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
    const char axisChar = axisName(axis)[0];
    DenseMatrix total(stateCount(kSpins));
    for (SpinIndex s = 0; s < kSpins; ++s) {
      total += pauliOperator(kSpins, s, axisChar);
    }
    const StateVector applied = denseApply(total, psi);
    const auto [wantRe, wantIm] = applied.innerProduct(ref);

    const TotalSpinCorrelator obs(axis);
    const std::vector<double> got = record(obs, contextFor(psi, ref, h));

    INFO("axis=" << axisName(axis));
    REQUIRE(got[0] == Catch::Approx(wantRe).margin(1e-13));
    REQUIRE(got[1] == Catch::Approx(wantIm).margin(1e-13));
  }
}

TEST_CASE("one correlator class covers every legacy directory variant",
          "[observable]") {
  // CPMGZ measured x against an x-prepared reference; Rabi_10_SzSz measured z
  // against a z-prepared one; Rabi_10_SxSz mixed the two. All three are the
  // same code here, differing only in an argument.
  constexpr std::size_t kSpins = 3;
  const StateVector psi = randomState(kSpins, 0x3333);
  const Hamiltonian h(kSpins);

  for (const Axis prepared : {Axis::X, Axis::Z}) {
    StateVector reference(kSpins);
    DenseMatrix total(stateCount(kSpins));
    for (SpinIndex s = 0; s < kSpins; ++s) {
      total += pauliOperator(kSpins, s, axisName(prepared)[0]);
    }
    reference = denseApply(total, psi);

    for (const Axis measured : {Axis::X, Axis::Z}) {
      const TotalSpinCorrelator obs(measured);
      const std::vector<double> got = record(obs, contextFor(psi, reference, h));
      INFO("prepared=" << axisName(prepared)
                       << " measured=" << axisName(measured));
      REQUIRE(std::isfinite(got[0]));
      REQUIRE(std::isfinite(got[1]));
    }
  }
}

TEST_CASE("magnetisation of a known state", "[observable]") {
  // All spins up: every sigma_z reads +1, sigma_x and sigma_y read 0.
  constexpr std::size_t kSpins = 3;
  StateVector psi(kSpins);
  psi.setBasisState(stateCount(kSpins) - 1);
  const Hamiltonian h(kSpins);
  MeasureContext ctx;
  ctx.state = &psi;
  ctx.hamiltonian = &h;

  for (SpinIndex s = 0; s < kSpins; ++s) {
    REQUIRE(record(SingleSpinMagnetization(s, Axis::Z), ctx)[0] ==
            Catch::Approx(1.0));
    REQUIRE(record(SingleSpinMagnetization(s, Axis::X), ctx)[0] ==
            Catch::Approx(0.0).margin(1e-15));
    REQUIRE(record(SingleSpinMagnetization(s, Axis::Y), ctx)[0] ==
            Catch::Approx(0.0).margin(1e-15));
  }
}

TEST_CASE("energy is the Hamiltonian expectation and is conserved",
          "[observable]") {
  constexpr std::size_t kSpins = 4;
  Hamiltonian h(kSpins);
  for (SpinIndex s = 0; s < kSpins; ++s) h.setField(s, Field{0.3, 0.1, 0.7});
  h.setCoupling(SpinPair{0, 1}, 0.5, 0.5, 1.0);
  h.setCoupling(SpinPair{1, 2}, 0.4, -0.2, 0.9);

  StateVector psi = randomState(kSpins, 0x4444);
  const StateVector ref = psi;

  const Energy obs;
  const double before = record(obs, contextFor(psi, ref, h))[0];

  Propagator prop(kSpins, 1e-12);
  prop.evolve(h, TimeMode::Real, psi, 0.05, 200);

  const double after = record(obs, contextFor(psi, ref, h))[0];
  REQUIRE(after == Catch::Approx(before).margin(1e-11));
}

TEST_CASE("norms report the propagated and reference states", "[observable]") {
  constexpr std::size_t kSpins = 3;
  StateVector psi = randomState(kSpins, 0x5555);
  StateVector ref = randomState(kSpins, 0x6666);
  ref.scale(2.0);
  const Hamiltonian h(kSpins);

  const std::vector<double> got =
      record(StateNorms(), contextFor(psi, ref, h));
  REQUIRE(got[0] == Catch::Approx(1.0));
  REQUIRE(got[1] == Catch::Approx(2.0));
}

TEST_CASE("observables describe their own output columns", "[observable]") {
  // The writer takes headings from here, so a results file says what it holds
  // -- unlike vecResOut, whose layout lived only in the reader's head.
  REQUIRE(TotalSpinCorrelator(Axis::Z).columns() ==
          std::vector<std::string>{"corr_z_re", "corr_z_im"});
  REQUIRE(SingleSpinMagnetization(0, Axis::X).columns() ==
          std::vector<std::string>{"mag_1_x"});
  REQUIRE(StateNorms().columns().size() == 2);
  REQUIRE(Energy().columns().size() == 1);

  REQUIRE(TotalSpinCorrelator(Axis::Y).width() == 2);
}

TEST_CASE("observables are built from configuration strings", "[observable]") {
  REQUIRE(makeObservable("corr:z")->name() == "corr:z");
  REQUIRE(makeObservable("corr:X")->name() == "corr:x");
  REQUIRE(makeObservable("norm")->name() == "norm");
  REQUIRE(makeObservable("energy")->name() == "energy");
  // Spins are 1-based in configuration, as in the legacy .dat files.
  REQUIRE(makeObservable("mag:2:y")->name() == "mag:2:y");
}

TEST_CASE("bad observable specifications are rejected with guidance",
          "[observable][errors]") {
  REQUIRE_THROWS_AS(makeObservable("corr:w"), std::invalid_argument);
  REQUIRE_THROWS_AS(makeObservable("corr"), std::invalid_argument);
  REQUIRE_THROWS_AS(makeObservable("mag:0:z"), std::invalid_argument);
  REQUIRE_THROWS_AS(makeObservable("mag:two:z"), std::invalid_argument);
  REQUIRE_THROWS_AS(makeObservable("nonsense"), std::invalid_argument);

  try {
    [[maybe_unused]] const auto ignored = makeObservable("nonsense");
  } catch (const std::invalid_argument& e) {
    const std::string message = e.what();
    REQUIRE(message.find("corr:<axis>") != std::string::npos);
  }
}

TEST_CASE("a correlator without a reference state is an error",
          "[observable][errors]") {
  const StateVector psi = randomState(2, 0x7777);
  MeasureContext ctx;
  ctx.state = &psi;

  std::vector<double> out(2);
  REQUIRE_THROWS_AS(TotalSpinCorrelator(Axis::Z).measure(ctx, out),
                    std::invalid_argument);
}
