#include "spinsim/core/StateVector.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace spinsim {

void StateVector::setZero() noexcept { std::ranges::fill(data_, 0.0); }

void StateVector::scale(double factor) noexcept {
  // The interleaved layout lets every one of these run over one flat array,
  // with no distinction between real and imaginary parts.
  for (double& v : data_) v *= factor;
}

void StateVector::negate() noexcept {
  for (double& v : data_) v = -v;
}

void StateVector::addScaled(double cr, double ci,
                            const StateVector& other) noexcept {
  const double* src = other.data_.data();
  double* dst = data_.data();
  const std::size_t n = data_.size();
  for (std::size_t k = 0; k < n; k += 2) {
    const double orr = src[k];
    const double oi = src[k + 1];
    dst[k] += cr * orr - ci * oi;
    dst[k + 1] += cr * oi + ci * orr;
  }
}

void StateVector::setBasisState(std::size_t index) {
  if (index >= size()) {
    throw std::out_of_range("basis state index outside the 2^L state space");
  }
  setZero();
  re(index) = 1.0;
}

double StateVector::norm() const noexcept {
  // Real and imaginary parts contribute identically to the sum of squares, so
  // this is one flat reduction.
  double sum = 0.0;
  for (const double v : data_) sum += v * v;
  return std::sqrt(sum);
}

double StateVector::normalize() {
  const double n = norm();
  if (n == 0.0) throw std::domain_error("cannot normalize a zero state");
  scale(1.0 / n);
  return n;
}

std::pair<double, double> StateVector::innerProduct(
    const StateVector& other) const noexcept {
  double real = 0.0;
  double imag = 0.0;
  const double* a = data_.data();
  const double* b = other.data_.data();
  const std::size_t n = data_.size();
  for (std::size_t k = 0; k < n; k += 2) {
    real += b[k] * a[k] + b[k + 1] * a[k + 1];
    imag += b[k] * a[k + 1] - b[k + 1] * a[k];
  }
  return {real, imag};
}

void StateVector::swap(StateVector& other) noexcept {
  std::swap(nspins_, other.nspins_);
  data_.swap(other.data_);
}

}  // namespace spinsim
