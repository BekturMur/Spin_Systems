#pragma once

#include <cstddef>
#include <cstdlib>
#include <limits>
#include <new>
#include <span>
#include <utility>
#include <vector>

#include "spinsim/core/Bits.hpp"

namespace spinsim {

/// Allocator handing out `Alignment`-aligned storage, so the amplitude array
/// starts on a cache-line boundary and the kernels vectorise cleanly.
template <typename T, std::size_t Alignment = 64>
struct AlignedAllocator {
  using value_type = T;

  AlignedAllocator() noexcept = default;
  template <typename U>
  explicit AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept {}

  template <typename U>
  struct rebind {
    using other = AlignedAllocator<U, Alignment>;
  };

  [[nodiscard]] T* allocate(std::size_t n) {
    if (n > std::numeric_limits<std::size_t>::max() / sizeof(T)) {
      throw std::bad_alloc{};
    }
    // std::aligned_alloc requires the size to be a multiple of the alignment.
    const std::size_t bytes =
        ((n * sizeof(T) + Alignment - 1) / Alignment) * Alignment;
    void* p = std::aligned_alloc(Alignment, bytes);
    if (p == nullptr) throw std::bad_alloc{};
    return static_cast<T*>(p);
  }

  void deallocate(T* p, std::size_t) noexcept { std::free(p); }

  template <typename U>
  friend bool operator==(const AlignedAllocator&,
                         const AlignedAllocator<U, Alignment>&) noexcept {
    return true;
  }
};

using AlignedDoubles = std::vector<double, AlignedAllocator<double>>;

/// A many-spin pure state: 2^L complex amplitudes, stored interleaved.
///
/// Amplitude k occupies data()[2k] and data()[2k + 1]. The legacy code kept two
/// separate arrays for the real and imaginary parts (`psiR`/`psiI`), and so did
/// this one until the layout was changed.
///
/// Be warned that the change bought nothing. In isolation the two-spin
/// flip-flop kernel -- four fifths of the run -- is 1.7x to 1.8x faster
/// interleaved, because it moves whole amplitudes between index pairs and the
/// split layout made that two unrelated streams. End to end the propagator did
/// not move at all: 0.412 ms per expansion order before, 0.414 ms after, at
/// fourteen spins.
///
/// The reason is that the isolated benchmark was not representative. It ran one
/// kernel against a single 256 KB input, which is compute-bound. The recurrence
/// keeps three state vectors and the diagonal live, close to a megabyte, and
/// there the run is limited by bytes moved rather than by instructions issued --
/// and both layouts move exactly the same bytes. Re-measuring the kernel with
/// comparable buffers live shrinks its advantage from 1.81x to 1.28x, and the
/// rest goes in the propagator's own passes.
///
/// The layout is kept because it is equivalent in speed, needs one allocation
/// instead of two, and reads better at the call sites. It is not kept because
/// it is faster. Anyone tempted to revisit this should benchmark the whole
/// propagator, not a kernel.
///
/// The length is the run's actual 2^L rather than the legacy `maxStat = 2^24`.
/// chstepsPDDGnmr.f declared six such arrays as locals, so a ten-spin run
/// reserved about 0.8 GB to use 16 KB of it.
class StateVector {
 public:
  StateVector() = default;

  /// Zero-initialised state over `nspins` spins.
  explicit StateVector(std::size_t nspins)
      : nspins_(nspins), data_(2 * stateCount(nspins), 0.0) {}

  [[nodiscard]] std::size_t nspins() const noexcept { return nspins_; }

  /// Number of complex amplitudes, i.e. 2^L. The backing array is twice this.
  [[nodiscard]] std::size_t size() const noexcept { return data_.size() / 2; }

  /// The interleaved backing store, of length 2 * size().
  [[nodiscard]] double* data() noexcept { return data_.data(); }
  [[nodiscard]] const double* data() const noexcept { return data_.data(); }
  [[nodiscard]] std::span<double> amplitudes() noexcept { return data_; }
  [[nodiscard]] std::span<const double> amplitudes() const noexcept {
    return data_;
  }

  [[nodiscard]] double& re(std::size_t k) noexcept { return data_[2 * k]; }
  [[nodiscard]] double& im(std::size_t k) noexcept { return data_[2 * k + 1]; }
  [[nodiscard]] double re(std::size_t k) const noexcept { return data_[2 * k]; }
  [[nodiscard]] double im(std::size_t k) const noexcept {
    return data_[2 * k + 1];
  }

  /// Set every amplitude to zero, keeping the allocation.
  void setZero() noexcept;

  /// Multiply every amplitude by a real factor.
  void scale(double factor) noexcept;

  /// Negate every amplitude, i.e. scale(-1) without the multiply.
  void negate() noexcept;

  /// this += (cr + i ci) * other. Both states must span the same space.
  void addScaled(double cr, double ci, const StateVector& other) noexcept;

  /// Collapse onto the single basis state `index`.
  void setBasisState(std::size_t index);

  /// Euclidean norm sqrt(sum |c_k|^2).
  [[nodiscard]] double norm() const noexcept;

  /// Rescale to unit norm. Returns the norm found before scaling.
  double normalize();

  /// Overlap <other|this> = sum conj(other_k) * this_k.
  [[nodiscard]] std::pair<double, double> innerProduct(
      const StateVector& other) const noexcept;

  void swap(StateVector& other) noexcept;

 private:
  std::size_t nspins_ = 0;
  AlignedDoubles data_;
};

}  // namespace spinsim
