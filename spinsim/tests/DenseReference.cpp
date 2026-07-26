#include "DenseReference.hpp"

#include <cmath>
#include <numbers>
#include <random>
#include <stdexcept>

namespace spinsim::test {

DenseMatrix DenseMatrix::identity(std::size_t dim) {
  DenseMatrix m(dim);
  for (std::size_t i = 0; i < dim; ++i) m(i, i) = Complex{1.0, 0.0};
  return m;
}

DenseMatrix& DenseMatrix::operator+=(const DenseMatrix& rhs) {
  if (dim_ != rhs.dim_) throw std::invalid_argument("dimension mismatch");
  for (std::size_t i = 0; i < data_.size(); ++i) data_[i] += rhs.data_[i];
  return *this;
}

DenseMatrix DenseMatrix::operator*(Complex scale) const {
  DenseMatrix out = *this;
  for (Complex& v : out.data_) v *= scale;
  return out;
}

DenseMatrix kron(const DenseMatrix& lhs, const DenseMatrix& rhs) {
  const std::size_t rd = rhs.dim();
  DenseMatrix out(lhs.dim() * rd);
  for (std::size_t lr = 0; lr < lhs.dim(); ++lr) {
    for (std::size_t lc = 0; lc < lhs.dim(); ++lc) {
      for (std::size_t rr = 0; rr < rd; ++rr) {
        for (std::size_t rc = 0; rc < rd; ++rc) {
          out(lr * rd + rr, lc * rd + rc) = lhs(lr, lc) * rhs(rr, rc);
        }
      }
    }
  }
  return out;
}

DenseMatrix pauli2x2(char axis) {
  DenseMatrix m(2);
  // Basis order is (|down>, |up>), i.e. index 0 is spin-down, matching the
  // convention that a set bit means "up".
  switch (axis) {
    case 'x':
      m(0, 1) = Complex{1.0, 0.0};
      m(1, 0) = Complex{1.0, 0.0};
      break;
    case 'y':
      // sigma_y |down> = -i|up>, sigma_y |up> = +i|down>. With index 0 = down
      // this is the usual [[0,-i],[i,0]] with its basis order reversed. The
      // sign is fixed by sigma_x sigma_y = i sigma_z, checked in the tests.
      m(0, 1) = Complex{0.0, 1.0};
      m(1, 0) = Complex{0.0, -1.0};
      break;
    case 'z':
      m(0, 0) = Complex{-1.0, 0.0};
      m(1, 1) = Complex{1.0, 0.0};
      break;
    default:
      throw std::invalid_argument("axis must be one of x, y, z");
  }
  return m;
}

DenseMatrix pauliOperator(std::size_t nspins, SpinIndex spin, char axis) {
  if (spin >= nspins) throw std::out_of_range("spin index out of range");
  // Spin s owns bit s, so it sits at position (nspins - 1 - s) counting from
  // the most significant Kronecker factor.
  DenseMatrix out = DenseMatrix::identity(1);
  for (std::size_t pos = nspins; pos-- > 0;) {
    out = kron(out, pos == spin ? pauli2x2(axis) : DenseMatrix::identity(2));
  }
  return out;
}

namespace {

/// Dense matrix product, only used to form sigma_a(i) * sigma_a(j).
DenseMatrix multiply(const DenseMatrix& a, const DenseMatrix& b) {
  const std::size_t n = a.dim();
  DenseMatrix out(n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t k = 0; k < n; ++k) {
      const Complex aik = a(i, k);
      if (aik == Complex{0.0, 0.0}) continue;
      for (std::size_t j = 0; j < n; ++j) out(i, j) += aik * b(k, j);
    }
  }
  return out;
}

}  // namespace

DenseMatrix couplingOperator(std::size_t nspins, SpinPair pair, double jx,
                             double jy, double jz) {
  DenseMatrix out(stateCount(nspins));
  const double j[3] = {jx, jy, jz};
  const char axes[3] = {'x', 'y', 'z'};
  for (int a = 0; a < 3; ++a) {
    if (j[static_cast<std::size_t>(a)] == 0.0) continue;
    const DenseMatrix term =
        multiply(pauliOperator(nspins, pair.low(), axes[a]),
                 pauliOperator(nspins, pair.high(), axes[a]));
    out += term * Complex{j[static_cast<std::size_t>(a)], 0.0};
  }
  return out;
}

DenseMatrix denseHamiltonian(const Hamiltonian& h) {
  const std::size_t n = h.nspins();
  DenseMatrix out(stateCount(n));
  // S = sigma/2 turns field terms into H/2 and coupling terms into J/4.
  for (SpinIndex s = 0; s < n; ++s) {
    const Field& f = h.field(s);
    const double c[3] = {f.x, f.y, f.z};
    const char axes[3] = {'x', 'y', 'z'};
    for (int a = 0; a < 3; ++a) {
      const double v = c[static_cast<std::size_t>(a)];
      if (v == 0.0) continue;
      out += pauliOperator(n, s, axes[a]) * Complex{0.5 * v, 0.0};
    }
  }
  for (const Coupling& c : h.couplings()) {
    out += couplingOperator(n, c.pair, 0.25 * c.jx, 0.25 * c.jy, 0.25 * c.jz);
  }
  return out;
}

StateVector denseApply(const DenseMatrix& m, const StateVector& in) {
  if (m.dim() != in.size()) throw std::invalid_argument("dimension mismatch");
  StateVector out(in.nspins());
  for (std::size_t r = 0; r < m.dim(); ++r) {
    Complex acc{0.0, 0.0};
    for (std::size_t c = 0; c < m.dim(); ++c) {
      acc += m(r, c) * Complex{in.re(c), in.im(c)};
    }
    out.re(r) = acc.real();
    out.im(r) = acc.imag();
  }
  return out;
}

extern "C" void zheev_(const char* jobz, const char* uplo, const int* n,
                       double* a, const int* lda, double* w, double* work,
                       const int* lwork, double* rwork, int* info);

StateVector denseExpApply(const DenseMatrix& m, Complex factor,
                          const StateVector& in) {
  const auto n = static_cast<int>(m.dim());

  // LAPACK wants column-major; copy across explicitly rather than relying on
  // the matrix being Hermitian to make the layouts coincide.
  std::vector<Complex> a(m.dim() * m.dim());
  for (std::size_t r = 0; r < m.dim(); ++r) {
    for (std::size_t c = 0; c < m.dim(); ++c) a[c * m.dim() + r] = m(r, c);
  }

  std::vector<double> eigenvalues(m.dim());
  std::vector<Complex> work(static_cast<std::size_t>(2 * n) * 2 + 2);
  std::vector<double> rwork(static_cast<std::size_t>(3 * n));
  const int lwork = static_cast<int>(work.size());
  int info = 0;
  zheev_("V", "U", &n, reinterpret_cast<double*>(a.data()), &n,
         eigenvalues.data(), reinterpret_cast<double*>(work.data()), &lwork,
         rwork.data(), &info);
  if (info != 0) throw std::runtime_error("zheev failed");

  // a now holds the eigenvectors as columns. Transform in -> eigenbasis, scale
  // by exp(factor * lambda), transform back.
  const std::size_t dim = m.dim();
  std::vector<Complex> coeff(dim, Complex{0.0, 0.0});
  for (std::size_t k = 0; k < dim; ++k) {
    Complex acc{0.0, 0.0};
    for (std::size_t r = 0; r < dim; ++r) {
      acc += std::conj(a[k * dim + r]) * Complex{in.re(r), in.im(r)};
    }
    coeff[k] = acc * std::exp(factor * eigenvalues[k]);
  }

  StateVector out(in.nspins());
  for (std::size_t r = 0; r < dim; ++r) {
    Complex acc{0.0, 0.0};
    for (std::size_t k = 0; k < dim; ++k) acc += a[k * dim + r] * coeff[k];
    out.re(r) = acc.real();
    out.im(r) = acc.imag();
  }
  return out;
}

double maxAbsDiff(const StateVector& a, const StateVector& b) {
  double worst = 0.0;
  for (std::size_t k = 0; k < a.size(); ++k) {
    const double dr = a.re(k) - b.re(k);
    const double di = a.im(k) - b.im(k);
    worst = std::max(worst, std::hypot(dr, di));
  }
  return worst;
}

StateVector randomState(std::size_t nspins, std::uint64_t seed) {
  StateVector psi(nspins);
  std::mt19937_64 rng(seed);
  std::normal_distribution<double> gauss(0.0, 1.0);
  for (std::size_t k = 0; k < psi.size(); ++k) {
    psi.re(k) = gauss(rng);
    psi.im(k) = gauss(rng);
  }
  psi.normalize();
  return psi;
}

}  // namespace spinsim::test
