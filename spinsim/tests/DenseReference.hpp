#pragma once

#include <complex>
#include <cstddef>
#include <vector>

#include "spinsim/core/Bits.hpp"
#include "spinsim/core/Hamiltonian.hpp"
#include "spinsim/core/StateVector.hpp"

/// Independent dense reference for the sparse kernels.
///
/// Everything here is built from 2x2 blocks via explicit Kronecker products and
/// dense matrix-vector products. Nothing reuses the bit arithmetic under test,
/// so agreement between the two is real evidence rather than a tautology.
/// Only usable for small L -- a dense 2^L x 2^L matrix is the whole point.
namespace spinsim::test {

using Complex = std::complex<double>;

/// Row-major dense square matrix.
class DenseMatrix {
 public:
  DenseMatrix() = default;
  explicit DenseMatrix(std::size_t dim)
      : dim_(dim), data_(dim * dim, Complex{0.0, 0.0}) {}

  [[nodiscard]] std::size_t dim() const noexcept { return dim_; }
  Complex& operator()(std::size_t r, std::size_t c) { return data_[r * dim_ + c]; }
  const Complex& operator()(std::size_t r, std::size_t c) const {
    return data_[r * dim_ + c];
  }

  static DenseMatrix identity(std::size_t dim);

  DenseMatrix& operator+=(const DenseMatrix& rhs);
  DenseMatrix operator*(Complex scale) const;

 private:
  std::size_t dim_ = 0;
  std::vector<Complex> data_;
};

/// Kronecker product, `lhs` occupying the more significant index bits.
[[nodiscard]] DenseMatrix kron(const DenseMatrix& lhs, const DenseMatrix& rhs);

/// The 2x2 Pauli matrices in the basis (index 0 = down, index 1 = up).
[[nodiscard]] DenseMatrix pauli2x2(char axis);

/// Single-spin Pauli operator embedded in the full 2^nspins space.
[[nodiscard]] DenseMatrix pauliOperator(std::size_t nspins, SpinIndex spin,
                                        char axis);

/// jx sx.sx + jy sy.sy + jz sz.sz for one pair, in the full space.
[[nodiscard]] DenseMatrix couplingOperator(std::size_t nspins, SpinPair pair,
                                           double jx, double jy, double jz);

/// Dense form of the physical (spin-1/2) Hamiltonian.
[[nodiscard]] DenseMatrix denseHamiltonian(const Hamiltonian& h);

/// out = m * in, using dense arithmetic.
[[nodiscard]] StateVector denseApply(const DenseMatrix& m,
                                     const StateVector& in);

/// exp(factor * m) * in, via eigendecomposition of the Hermitian `m`.
///
/// This is the yardstick the Chebyshev propagator is measured against: an
/// entirely different algorithm for the same quantity, exact up to LAPACK's own
/// accuracy.
[[nodiscard]] StateVector denseExpApply(const DenseMatrix& m, Complex factor,
                                        const StateVector& in);

/// Largest |a_k - b_k| over all amplitudes.
[[nodiscard]] double maxAbsDiff(const StateVector& a, const StateVector& b);

/// Deterministic pseudo-random normalised state, for exercising every amplitude.
[[nodiscard]] StateVector randomState(std::size_t nspins, std::uint64_t seed);

}  // namespace spinsim::test
