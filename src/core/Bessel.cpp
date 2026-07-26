#include "spinsim/core/Bessel.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace spinsim {
namespace {

/// Rescaling thresholds for the backward recurrence. Terms grow rapidly as the
/// index falls below the argument, so the whole partial result is scaled down
/// whenever it approaches the top of the double range.
constexpr double kBig = 1.0e250;
constexpr double kInvBig = 1.0e-250;

/// Where to start the downward recurrence. It must begin far enough above both
/// the requested order and the argument that the unwanted second solution has
/// decayed away by the time the recurrence reaches index 0. The square-root
/// margin is the standard choice; the additive term keeps small cases sane.
[[nodiscard]] std::size_t startIndex(double x, std::size_t count) {
  const double top = std::max(static_cast<double>(count), std::ceil(x));
  const double margin = 40.0 + 6.0 * std::sqrt(top);
  const auto start = static_cast<std::size_t>(top + margin);
  return start + (start % 2);  // even, so the J normalisation sum ends cleanly
}

/// Trims trailing zeros and records how many leading entries survived.
[[nodiscard]] std::size_t countNonzero(const std::vector<double>& v) {
  std::size_t n = v.size();
  while (n > 0 && v[n - 1] == 0.0) --n;
  return n;
}

void validate(double x, std::size_t count) {
  if (count == 0) {
    throw std::invalid_argument("Bessel sequence length must be positive");
  }
  if (!(x >= 0.0)) {
    throw std::invalid_argument("Bessel argument must be non-negative");
  }
}

}  // namespace

BesselSequence besselJ(double x, std::size_t count) {
  validate(x, count);

  BesselSequence out;
  out.values.assign(count, 0.0);

  if (x == 0.0) {
    out.values[0] = 1.0;
    out.nonzero = 1;
    return out;
  }

  // Backward recurrence  J_{k-1} = (2k/x) J_k - J_{k+1}, seeded with an
  // arbitrary tiny value; the result is correct only up to a common factor,
  // which the normalisation below removes.
  const std::size_t start = startIndex(x, count);
  double higher = 0.0;      // stands for J_{k+1}
  double current = 1.0e-30; // stands for J_k
  double normalization = 0.0;

  for (std::size_t k = start; k >= 1; --k) {
    const double lower = (2.0 * static_cast<double>(k) / x) * current - higher;
    higher = current;
    current = lower;

    // Every even order except zero enters the identity twice.
    if (k % 2 == 1) normalization += 2.0 * current;

    if (std::abs(current) > kBig) {
      current *= kInvBig;
      higher *= kInvBig;
      normalization *= kInvBig;
      for (std::size_t i = 0; i < count; ++i) out.values[i] *= kInvBig;
    }

    // `current` now holds order k-1; keep it if it is one we were asked for.
    if (k - 1 < count) out.values[k - 1] = current;
  }

  // J_0 + 2(J_2 + J_4 + ...) = 1. The loop above accumulated the even orders
  // (index k-1 even, i.e. k odd) including order zero, so subtract it once.
  normalization -= current;

  if (normalization == 0.0 || !std::isfinite(normalization)) {
    throw std::runtime_error("Bessel J recurrence failed to normalise");
  }
  const double scale = 1.0 / normalization;
  for (double& v : out.values) v *= scale;

  out.nonzero = countNonzero(out.values);
  return out;
}

BesselSequence besselI(double x, std::size_t count) {
  validate(x, count);

  BesselSequence out;
  out.values.assign(count, 0.0);

  if (x == 0.0) {
    out.values[0] = 1.0;
    out.nonzero = 1;
    return out;
  }
  // I_0(x) alone already exceeds the double range beyond here.
  if (x > 709.0) {
    throw std::overflow_error("unscaled modified Bessel I overflows");
  }

  // Same scheme, with the sign in the recurrence flipped:
  //   I_{k-1} = (2k/x) I_k + I_{k+1}
  const std::size_t start = startIndex(x, count);
  double higher = 0.0;
  double current = 1.0e-30;
  double normalization = 0.0;

  for (std::size_t k = start; k >= 1; --k) {
    const double lower = (2.0 * static_cast<double>(k) / x) * current + higher;
    higher = current;
    current = lower;

    normalization += 2.0 * current;  // every order but zero, corrected below

    if (std::abs(current) > kBig) {
      current *= kInvBig;
      higher *= kInvBig;
      normalization *= kInvBig;
      for (std::size_t i = 0; i < count; ++i) out.values[i] *= kInvBig;
    }

    if (k - 1 < count) out.values[k - 1] = current;
  }

  // I_0 + 2(I_1 + I_2 + ...) = exp(x); order zero was counted twice above.
  normalization -= current;

  if (normalization == 0.0 || !std::isfinite(normalization)) {
    throw std::runtime_error("Bessel I recurrence failed to normalise");
  }
  const double scale = std::exp(x) / normalization;
  for (double& v : out.values) v *= scale;

  out.nonzero = countNonzero(out.values);
  return out;
}

}  // namespace spinsim
