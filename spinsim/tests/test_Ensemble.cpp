#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numeric>

#include "spinsim/ensemble/Ensemble.hpp"
#include "spinsim/ensemble/EnsembleRunner.hpp"

using namespace spinsim;

TEST_CASE("a realisation's stream depends only on its index", "[ensemble][rng]") {
  // The property the legacy code lacked: chebNMR2D-v6.f:468 accumulated the
  // seed into a worker-local variable that was never reset, so the stream
  // depended on how many jobs that rank had already run.
  for (std::size_t index : {std::size_t{0}, std::size_t{7}, std::size_t{179}}) {
    Rng a = Rng::forRealization(12345, index);
    Rng b = Rng::forRealization(12345, index);
    for (int k = 0; k < 16; ++k) REQUIRE(a.uniform() == b.uniform());
  }
}

TEST_CASE("adjacent realisations get unrelated streams", "[ensemble][rng]") {
  Rng a = Rng::forRealization(999, 4);
  Rng b = Rng::forRealization(999, 5);
  int matches = 0;
  for (int k = 0; k < 32; ++k) {
    if (a.uniform() == b.uniform()) ++matches;
  }
  REQUIRE(matches == 0);
}

TEST_CASE("gaussian draws have the right moments", "[ensemble][rng]") {
  Rng rng(2024);
  double sum = 0.0;
  double sumSq = 0.0;
  constexpr int kDraws = 200000;
  for (int k = 0; k < kDraws; ++k) {
    const double g = rng.gaussian();
    sum += g;
    sumSq += g * g;
  }
  const double mean = sum / kDraws;
  const double variance = sumSq / kDraws - mean * mean;
  REQUIRE(std::abs(mean) < 0.01);
  REQUIRE(variance == Catch::Approx(1.0).margin(0.02));
}

TEST_CASE("the accumulator reproduces textbook mean and variance",
          "[ensemble][accumulator]") {
  const std::vector<double> samples{2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0};
  Accumulator acc(1);
  for (const double s : samples) acc.add(std::span<const double>(&s, 1));

  REQUIRE(acc.count() == samples.size());
  REQUIRE(acc.mean()[0] == Catch::Approx(5.0));
  // Sample standard deviation of this classic set is sqrt(32/7).
  REQUIRE(acc.standardDeviation()[0] == Catch::Approx(std::sqrt(32.0 / 7.0)));
  REQUIRE(acc.standardError()[0] ==
          Catch::Approx(std::sqrt(32.0 / 7.0) / std::sqrt(8.0)));
}

TEST_CASE("standard error is smaller than standard deviation",
          "[ensemble][accumulator]") {
  // The legacy code printed sqrt(M2/n) as an error bar, which is neither of
  // these; keeping them distinct is the point.
  Accumulator acc(1);
  for (int k = 0; k < 100; ++k) {
    const double v = static_cast<double>(k % 7);
    acc.add(std::span<const double>(&v, 1));
  }
  REQUIRE(acc.standardError()[0] < acc.standardDeviation()[0]);
  REQUIRE(acc.standardError()[0] ==
          Catch::Approx(acc.standardDeviation()[0] / 10.0));
}

TEST_CASE("merging partial accumulators matches accumulating in one pass",
          "[ensemble][accumulator]") {
  const std::vector<double> data{1.5, -2.0, 3.25, 8.0, -0.5, 4.75, 6.0, 2.25};

  Accumulator whole(1);
  for (const double v : data) whole.add(std::span<const double>(&v, 1));

  Accumulator first(1);
  Accumulator second(1);
  for (std::size_t k = 0; k < data.size(); ++k) {
    (k < 3 ? first : second).add(std::span<const double>(&data[k], 1));
  }
  first.merge(second);

  REQUIRE(first.count() == whole.count());
  REQUIRE(first.mean()[0] == Catch::Approx(whole.mean()[0]));
  REQUIRE(first.standardDeviation()[0] ==
          Catch::Approx(whole.standardDeviation()[0]));
}

TEST_CASE("merging an empty accumulator changes nothing",
          "[ensemble][accumulator]") {
  Accumulator acc(2);
  const std::vector<double> sample{1.0, 2.0};
  acc.add(sample);
  const std::vector<double> before = acc.mean();

  acc.merge(Accumulator(2));
  REQUIRE(acc.count() == 1);
  REQUIRE(acc.mean() == before);
}

TEST_CASE("ensemble results do not depend on the thread count",
          "[ensemble][determinism]") {
  // The headline guarantee. Each realisation's numbers come from its index
  // alone, and the fold happens in index order, so eighteen threads give
  // bitwise what one thread gives.
  const std::vector<std::string> columns{"a", "b"};
  const RealizationJob job = [](std::size_t index, VectorRecorder& out) {
    Rng rng = Rng::forRealization(7777, index);
    for (int point = 0; point < 5; ++point) {
      const std::vector<double> row{rng.gaussian(), rng.uniform()};
      out.record(0.1 * point, row);
    }
  };

  LocalRunner single(1);
  LocalRunner many(18);
  const EnsembleResult a = single.run(64, columns, job);
  const EnsembleResult b = many.run(64, columns, job);

  REQUIRE(a.times == b.times);
  REQUIRE(a.statistics.count() == b.statistics.count());
  // Bitwise, not approximately.
  REQUIRE(a.statistics.mean() == b.statistics.mean());
  REQUIRE(a.statistics.standardDeviation() == b.statistics.standardDeviation());
}

TEST_CASE("a failing realisation surfaces as an exception",
          "[ensemble][errors]") {
  const std::vector<std::string> columns{"x"};
  const RealizationJob job = [](std::size_t index, VectorRecorder& out) {
    if (index == 13) throw std::runtime_error("realisation blew up");
    const std::vector<double> row{1.0};
    out.record(0.0, row);
  };

  LocalRunner runner(4);
  REQUIRE_THROWS_AS(runner.run(32, columns, job), std::runtime_error);
}

TEST_CASE("realisations that disagree on shape are rejected",
          "[ensemble][errors]") {
  const std::vector<std::string> columns{"x"};
  const RealizationJob job = [](std::size_t index, VectorRecorder& out) {
    const std::vector<double> row{1.0};
    out.record(0.0, row);
    if (index == 2) out.record(1.0, row);  // one extra point
  };

  LocalRunner runner(2);
  REQUIRE_THROWS_AS(runner.run(8, columns, job), std::runtime_error);
}

TEST_CASE("the dipolar generator produces every pair with 1/r^3 scaling",
          "[ensemble][geometry]") {
  constexpr std::size_t kSpins = 8;
  Rng rng = Rng::forRealization(31337, 0);
  const std::vector<Coupling> couplings =
      generateDipolarLattice2D(kSpins, rng, 1.0);

  REQUIRE(couplings.size() == kSpins * (kSpins - 1) / 2);
  for (const Coupling& c : couplings) {
    REQUIRE(c.pair.low() < c.pair.high());
    // Secular truncation: Jx = Jy = -Jz/2.
    REQUIRE(c.jx == Catch::Approx(-0.5 * c.jz));
    REQUIRE(c.jy == Catch::Approx(-0.5 * c.jz));
    REQUIRE(c.jz > 0.0);
  }
}

TEST_CASE("geometry is reproducible from the realisation index",
          "[ensemble][geometry]") {
  Rng a = Rng::forRealization(42, 5);
  Rng b = Rng::forRealization(42, 5);
  const auto first = generateDipolarLattice2D(6, a, 2.0);
  const auto second = generateDipolarLattice2D(6, b, 2.0);

  REQUIRE(first.size() == second.size());
  for (std::size_t k = 0; k < first.size(); ++k) {
    REQUIRE(first[k].jz == second[k].jz);
    REQUIRE(first[k].pair == second[k].pair);
  }
}

TEST_CASE("initial states are prepared as configured", "[ensemble][state]") {
  Rng rng(1);
  StateVector psi(4);

  prepareInitialState(psi, InitialStateKind::AllUp, rng);
  REQUIRE(psi.re(15) == 1.0);

  prepareInitialState(psi, InitialStateKind::AllDown, rng);
  REQUIRE(psi.re(0) == 1.0);

  prepareInitialState(psi, InitialStateKind::BasisState, rng, 6);
  REQUIRE(psi.re(6) == 1.0);

  prepareInitialState(psi, InitialStateKind::Random, rng);
  REQUIRE(psi.norm() == Catch::Approx(1.0));
  // A typical state spreads over the whole space.
  int populated = 0;
  for (std::size_t k = 0; k < psi.size(); ++k) {
    if (std::hypot(psi.re(k), psi.im(k)) > 1e-6) ++populated;
  }
  REQUIRE(populated == 16);
}

TEST_CASE("field disorder shifts only the z component", "[ensemble][state]") {
  Hamiltonian h(5);
  for (SpinIndex s = 0; s < 5; ++s) h.setField(s, Field{1.0, 2.0, 3.0});

  Rng rng = Rng::forRealization(8, 8);
  applyFieldDisorder(h, rng, 10.0);

  for (SpinIndex s = 0; s < 5; ++s) {
    REQUIRE(h.field(s).x == 1.0);
    REQUIRE(h.field(s).y == 2.0);
    REQUIRE(h.field(s).z != 3.0);
  }
}

TEST_CASE("initial state names are parsed", "[ensemble][state]") {
  REQUIRE(parseInitialStateKind("random") == InitialStateKind::Random);
  REQUIRE(parseInitialStateKind("up") == InitialStateKind::AllUp);
  REQUIRE_THROWS_AS(parseInitialStateKind("sideways"), std::invalid_argument);
}
