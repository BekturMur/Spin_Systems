#include "spinsim/core/Gate.hpp"

#include <charconv>
#include <cmath>
#include <format>
#include <numbers>
#include <stdexcept>

namespace spinsim {
namespace {

using Complex = std::complex<double>;

/// The 2x2 matrix exp(-i (theta/2) n.sigma) in the basis (down, up).
///
/// With sigma_z = diag(-1, +1) in that order, n.sigma is
///     [[-nz, nx + i ny], [nx - i ny, nz]],
/// and U = cos(theta/2) I - i sin(theta/2) (n.sigma). Written out this
/// reproduces the arithmetic of chspecstepsPDDGnmr.f:243-246 exactly.
struct Rotation2x2 {
  Complex ll, lh, hl, hh;

  static Rotation2x2 make(Axis3 n, double angle) {
    const double c = std::cos(0.5 * angle);
    const double s = std::sin(0.5 * angle);
    return Rotation2x2{
        Complex{c, s * n.z},
        Complex{s * n.y, -s * n.x},
        Complex{-s * n.y, -s * n.x},
        Complex{c, -s * n.z},
    };
  }
};

}  // namespace

Axis3 Axis3::normalized() const {
  const double length = std::sqrt(x * x + y * y + z * z);
  if (!(length > 0.0)) {
    throw std::invalid_argument("rotation axis must have non-zero length");
  }
  return Axis3{x / length, y / length, z / length};
}

SpinRotation::SpinRotation(Axis3 axis, double angle,
                           std::vector<SpinIndex> targets)
    : axis_(axis.normalized()), angle_(angle), targets_(std::move(targets)) {
  if (targets_.empty()) {
    throw std::invalid_argument("a rotation must name at least one spin");
  }
}

SpinRotation SpinRotation::onAllSpins(Axis3 axis, double angle,
                                      std::size_t nspins) {
  std::vector<SpinIndex> all(nspins);
  for (std::size_t s = 0; s < nspins; ++s) all[s] = s;
  return SpinRotation(axis, angle, std::move(all));
}

void SpinRotation::apply(StateVector& psi) const {
  const Rotation2x2 u = Rotation2x2::make(axis_, angle_);
  const std::size_t nstates = psi.size();
  double* d = psi.data();

  for (const SpinIndex spin : targets_) {
    if (spin >= psi.nspins()) {
      throw std::out_of_range("rotation names a spin outside the system");
    }
    const std::size_t stride = spinMask(spin);
    const std::size_t block = stride << 1;
    for (std::size_t base = 0; base < nstates; base += block) {
      for (std::size_t lo = base; lo < base + stride; ++lo) {
        const std::size_t hi = lo + stride;
        const Complex a{d[2 * lo], d[2 * lo + 1]};
        const Complex b{d[2 * hi], d[2 * hi + 1]};
        const Complex newLo = u.ll * a + u.lh * b;
        const Complex newHi = u.hl * a + u.hh * b;
        d[2 * lo] = newLo.real();
        d[2 * lo + 1] = newLo.imag();
        d[2 * hi] = newHi.real();
        d[2 * hi + 1] = newHi.imag();
      }
    }
  }
}

TwoSpinGate::TwoSpinGate(Matrix matrix, SpinPair pair)
    : matrix_(matrix), pair_(pair) {}

void TwoSpinGate::apply(StateVector& psi) const {
  if (pair_.high() >= psi.nspins()) {
    throw std::out_of_range("gate names a spin outside the system");
  }
  const std::size_t hiMask = spinMask(pair_.high());
  const std::size_t loMask = spinMask(pair_.low());
  const std::size_t hiBlock = hiMask << 1;
  const std::size_t loBlock = loMask << 1;
  const std::size_t nstates = psi.size();

  double* d = psi.data();

  for (std::size_t outer = 0; outer < nstates; outer += hiBlock) {
    for (std::size_t mid = outer; mid < outer + hiMask; mid += loBlock) {
      for (std::size_t dd = mid; dd < mid + loMask; ++dd) {
        // Local basis order (dd, du, ud, uu) with du = low spin up.
        const std::array<std::size_t, 4> idx{dd, dd + loMask, dd + hiMask,
                                             dd + hiMask + loMask};
        std::array<Complex, 4> in{};
        for (std::size_t k = 0; k < 4; ++k) {
          in[k] = Complex{d[2 * idx[k]], d[2 * idx[k] + 1]};
        }
        for (std::size_t r = 0; r < 4; ++r) {
          Complex acc{0.0, 0.0};
          for (std::size_t c = 0; c < 4; ++c) acc += matrix_[r * 4 + c] * in[c];
          d[2 * idx[r]] = acc.real();
          d[2 * idx[r] + 1] = acc.imag();
        }
      }
    }
  }
}

double TwoSpinGate::unitarityDefect() const {
  double worst = 0.0;
  for (std::size_t r = 0; r < 4; ++r) {
    for (std::size_t c = 0; c < 4; ++c) {
      Complex acc{0.0, 0.0};
      for (std::size_t k = 0; k < 4; ++k) {
        acc += std::conj(matrix_[k * 4 + r]) * matrix_[k * 4 + c];
      }
      const Complex expected = (r == c) ? Complex{1.0, 0.0} : Complex{0.0, 0.0};
      worst = std::max(worst, std::abs(acc - expected));
    }
  }
  return worst;
}

namespace {

std::vector<std::string_view> split(std::string_view spec) {
  std::vector<std::string_view> parts;
  std::size_t start = 0;
  while (true) {
    const std::size_t sep = spec.find(':', start);
    if (sep == std::string_view::npos) {
      parts.push_back(spec.substr(start));
      break;
    }
    parts.push_back(spec.substr(start, sep - start));
    start = sep + 1;
  }
  return parts;
}

double parseAngle(std::string_view text) {
  if (text == "pi") return std::numbers::pi;
  if (text == "pi/2") return 0.5 * std::numbers::pi;
  if (text == "pi/4") return 0.25 * std::numbers::pi;
  throw std::invalid_argument(std::format(
      "unknown pulse angle '{}'; expected pi, pi/2 or pi/4", text));
}

Axis3 parseRotationAxis(std::string_view text) {
  if (text == "x") return Axis3{1.0, 0.0, 0.0};
  if (text == "y") return Axis3{0.0, 1.0, 0.0};
  if (text == "z") return Axis3{0.0, 0.0, 1.0};
  throw std::invalid_argument(
      std::format("unknown rotation axis '{}'; expected x, y or z", text));
}

}  // namespace

SpinRotation makeRotation(std::string_view spec, std::size_t nspins) {
  const std::vector<std::string_view> parts = split(spec);
  if (parts.size() != 2 && parts.size() != 3) {
    throw std::invalid_argument(std::format(
        "malformed pulse '{}'; expected <angle>:<axis> or <angle>:<axis>:<spin>",
        spec));
  }

  const double angle = parseAngle(parts[0]);
  const Axis3 axis = parseRotationAxis(parts[1]);

  if (parts.size() == 2) return SpinRotation::onAllSpins(axis, angle, nspins);

  unsigned long index = 0;
  const auto* first = parts[2].data();
  const auto* last = first + parts[2].size();
  const auto [ptr, ec] = std::from_chars(first, last, index);
  if (ec != std::errc{} || ptr != last || index == 0 || index > nspins) {
    throw std::invalid_argument(std::format(
        "pulse '{}': spin must be between 1 and {}", spec, nspins));
  }
  return SpinRotation(axis, angle, {static_cast<SpinIndex>(index - 1)});
}

}  // namespace spinsim
