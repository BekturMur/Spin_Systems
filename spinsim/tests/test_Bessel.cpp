#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <string_view>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "spinsim/core/Bessel.hpp"

using namespace spinsim;

namespace {

/// One record of the reference file produced by gen_bessel_reference.f, which
/// calls the legacy Cody routines directly.
struct Reference {
  char kind = 'J';  // 'J' for rjbesl, 'I' for ribesl
  double x = 0.0;
  std::vector<double> values;
};

/// Parses one reference value, insisting the whole line is consumed.
///
/// std::stod is unusable here: it throws out_of_range on underflow, and the
/// deep tail of these sequences is full of subnormals and true zeros. strtod
/// reports the same condition through errno while still returning the right
/// value. Requiring full consumption catches malformed output -- notably the
/// Fortran habit of dropping the 'E' from three-digit exponents, which a
/// lenient parser reads as a number a hundred orders of magnitude too large.
double parseStrict(const std::string& line) {
  errno = 0;
  const char* begin = line.c_str();
  char* end = nullptr;
  const double value = std::strtod(begin, &end);

  REQUIRE(end != begin);
  REQUIRE(std::string_view(end).find_first_not_of(" \t\r") ==
          std::string_view::npos);
  // ERANGE on underflow is expected and harmless; on overflow it is not.
  REQUIRE(std::abs(value) < 1.0e300);
  return value;
}

std::vector<Reference> loadReference() {
  const std::string path = std::string(SPINSIM_TEST_DATA_DIR) +
                           "/bessel_reference.txt";
  std::ifstream in(path);
  REQUIRE(in.good());

  std::vector<Reference> records;
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream header(line);
    std::string kind;
    double x = 0.0;
    long nb = 0;
    long ncalc = 0;
    if (!(header >> kind >> x >> nb >> ncalc)) continue;

    Reference rec;
    rec.kind = kind[0];
    rec.x = x;
    rec.values.reserve(static_cast<std::size_t>(nb));
    for (long i = 0; i < nb; ++i) {
      REQUIRE(std::getline(in, line));
      INFO("reference line: " << line);
      rec.values.push_back(parseStrict(line));
    }
    records.push_back(std::move(rec));
  }
  return records;
}

/// Expansion coefficients matter relative to the size of the whole sequence,
/// so errors are measured against its largest entry rather than term by term.
double sequenceScale(const std::vector<double>& v) {
  double s = 0.0;
  for (double e : v) s = std::max(s, std::abs(e));
  return s;
}

}  // namespace

TEST_CASE("besselJ reproduces the legacy rjbesl.f", "[bessel]") {
  const auto records = loadReference();
  REQUIRE_FALSE(records.empty());

  std::size_t checked = 0;
  for (const Reference& rec : records) {
    if (rec.kind != 'J') continue;
    ++checked;

    const BesselSequence got = besselJ(rec.x, rec.values.size());
    REQUIRE(got.values.size() == rec.values.size());

    const double scale = sequenceScale(rec.values);
    for (std::size_t n = 0; n < rec.values.size(); ++n) {
      INFO("x=" << rec.x << " order=" << n << " got=" << got.values[n]
                << " want=" << rec.values[n]);
      REQUIRE(std::abs(got.values[n] - rec.values[n]) < 1e-14 * scale);
    }
  }
  REQUIRE(checked > 0);
}

TEST_CASE("besselI reproduces the legacy ribesl.f", "[bessel]") {
  const auto records = loadReference();

  std::size_t checked = 0;
  for (const Reference& rec : records) {
    if (rec.kind != 'I') continue;
    // ribesl was asked for unscaled values, which overflow past this point;
    // besselI reports that rather than returning infinities.
    if (rec.x > 709.0) continue;
    ++checked;

    const BesselSequence got = besselI(rec.x, rec.values.size());
    const double scale = sequenceScale(rec.values);
    for (std::size_t n = 0; n < rec.values.size(); ++n) {
      INFO("x=" << rec.x << " order=" << n << " got=" << got.values[n]
                << " want=" << rec.values[n]);
      REQUIRE(std::abs(got.values[n] - rec.values[n]) < 1e-13 * scale);
    }
  }
  REQUIRE(checked > 0);
}

TEST_CASE("Bessel functions satisfy their defining recurrences", "[bessel]") {
  // Independent of the reference data: J_{n-1} + J_{n+1} = (2n/x) J_n and
  // I_{n-1} - I_{n+1} = (2n/x) I_n.
  for (const double x : {0.25, 1.0, 4.0, 17.0, 60.0}) {
    const BesselSequence j = besselJ(x, 80);
    for (std::size_t n = 1; n + 1 < j.nonzero; ++n) {
      const double lhs = j.values[n - 1] + j.values[n + 1];
      const double rhs = (2.0 * static_cast<double>(n) / x) * j.values[n];
      INFO("J x=" << x << " n=" << n);
      REQUIRE(std::abs(lhs - rhs) < 1e-12 * std::max(1.0, std::abs(rhs)));
    }

    const BesselSequence i = besselI(x, 80);
    const double iScale = sequenceScale(i.values);
    for (std::size_t n = 1; n + 1 < i.nonzero; ++n) {
      const double lhs = i.values[n - 1] - i.values[n + 1];
      const double rhs = (2.0 * static_cast<double>(n) / x) * i.values[n];
      INFO("I x=" << x << " n=" << n);
      REQUIRE(std::abs(lhs - rhs) < 1e-12 * iScale);
    }
  }
}

TEST_CASE("Bessel normalisation identities hold", "[bessel]") {
  for (const double x : {0.5, 3.0, 12.0, 45.0}) {
    const BesselSequence j = besselJ(x, 300);
    double sum = j.values[0];
    for (std::size_t n = 2; n < j.values.size(); n += 2) sum += 2.0 * j.values[n];
    INFO("J identity at x=" << x);
    REQUIRE(sum == Catch::Approx(1.0).margin(1e-13));

    const BesselSequence i = besselI(x, 300);
    double isum = i.values[0];
    for (std::size_t n = 1; n < i.values.size(); ++n) isum += 2.0 * i.values[n];
    INFO("I identity at x=" << x);
    REQUIRE(isum == Catch::Approx(std::exp(x)).epsilon(1e-13));
  }
}

TEST_CASE("Bessel at zero argument", "[bessel]") {
  const BesselSequence j = besselJ(0.0, 5);
  REQUIRE(j.values[0] == 1.0);
  REQUIRE(j.nonzero == 1);
  for (std::size_t n = 1; n < 5; ++n) REQUIRE(j.values[n] == 0.0);

  const BesselSequence i = besselI(0.0, 5);
  REQUIRE(i.values[0] == 1.0);
  REQUIRE(i.nonzero == 1);
}

TEST_CASE("nonzero marks where the sequence underflows", "[bessel]") {
  // At tiny argument the coefficients fall off like (x/2)^n / n!, so only the
  // first few dozen orders are representable. The legacy routine signalled the
  // same thing through NCALC.
  const BesselSequence j = besselJ(1.0e-6, 400);
  REQUIRE(j.nonzero > 0);
  REQUIRE(j.nonzero < 400);
  for (std::size_t n = j.nonzero; n < 400; ++n) REQUIRE(j.values[n] == 0.0);
}

TEST_CASE("Bessel rejects invalid arguments", "[bessel]") {
  REQUIRE_THROWS_AS(besselJ(-1.0, 4), std::invalid_argument);
  REQUIRE_THROWS_AS(besselJ(1.0, 0), std::invalid_argument);
  REQUIRE_THROWS_AS(besselI(1000.0, 4), std::overflow_error);
}
