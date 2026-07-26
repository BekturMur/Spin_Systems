#include "spinsim/ensemble/Rng.hpp"

#include <cmath>
#include <numbers>

namespace spinsim {

double Rng::gaussian() {
  if (hasSpare_) {
    hasSpare_ = false;
    return spare_;
  }
  // Guard the logarithm against an exact zero, as chinitPDDGnmr.f:67 does with
  // its MDLT offset.
  constexpr double kFloor = 1.0e-15;
  const double u = kFloor + (1.0 - 2.0 * kFloor) * uniform();
  const double v = uniform();
  const double radius = std::sqrt(-2.0 * std::log(u));
  const double phase = 2.0 * std::numbers::pi * v;
  spare_ = radius * std::sin(phase);
  hasSpare_ = true;
  return radius * std::cos(phase);
}

}  // namespace spinsim
