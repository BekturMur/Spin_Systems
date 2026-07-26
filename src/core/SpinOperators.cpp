#include "spinsim/core/SpinOperators.hpp"

#include <cassert>

namespace spinsim {
namespace {

/// Checks shared by every kernel: same state space, and no aliasing.
///
/// These are preconditions the callers inside this library always satisfy, so
/// they are asserts rather than throws and vanish in release builds -- hence
/// [[maybe_unused]] on every parameter.
void checkOperands([[maybe_unused]] const StateVector& in,
                   [[maybe_unused]] const StateVector& out,
                   [[maybe_unused]] SpinIndex spin) {
  assert(in.size() == out.size() && "state vectors span different spaces");
  assert(&in != &out && "kernels require distinct input and output");
  assert(spin < in.nspins() && "spin index out of range");
}

void checkPair([[maybe_unused]] const StateVector& in,
               [[maybe_unused]] const StateVector& out,
               [[maybe_unused]] SpinPair pair) {
  assert(in.size() == out.size() && "state vectors span different spaces");
  assert(&in != &out && "kernels require distinct input and output");
  assert(pair.high() < in.nspins() && "spin index out of range");
}

/// Invokes `body(lo, hi)` for every pair of basis indices that differ only in
/// bit `spin`, with `lo` the "down" partner and `hi` the "up" partner.
///
/// The indices are amplitude indices; kernels double them to reach the
/// interleaved storage. This is the loop nest the legacy kernels each spelled
/// out by hand, with `ip0 = 2**(i0-1)` recomputed via integer exponentiation
/// inside the loop.
template <typename Body>
void forEachSpinPairIndex(std::size_t nstates, SpinIndex spin, Body body) {
  const std::size_t stride = spinMask(spin);
  const std::size_t block = stride << 1;
  for (std::size_t base = 0; base < nstates; base += block) {
    for (std::size_t lo = base; lo < base + stride; ++lo) {
      body(lo, lo + stride);
    }
  }
}

/// Invokes `body(first, second)` over the two corners a two-spin channel
/// connects, chosen by their offsets from the both-down index.
///
/// Kept separate from the four-corner walk below: a channel kernel that only
/// needs two corners should not pay for computing the other two. Measured at
/// 395 us against 520 us for ninety-one pairs at fourteen spins, which is 1.32x
/// -- enough to matter, since these kernels are four fifths of the run.
template <typename Body>
void forEachCouplingChannel(std::size_t nstates, SpinPair pair,
                            std::size_t firstOffset, std::size_t secondOffset,
                            Body body) {
  const std::size_t hiMask = spinMask(pair.high());
  const std::size_t loMask = spinMask(pair.low());
  const std::size_t hiBlock = hiMask << 1;
  const std::size_t loBlock = loMask << 1;
  for (std::size_t outer = 0; outer < nstates; outer += hiBlock) {
    for (std::size_t mid = outer; mid < outer + hiMask; mid += loBlock) {
      for (std::size_t dd = mid; dd < mid + loMask; ++dd) {
        body(2 * (dd + firstOffset), 2 * (dd + secondOffset));
      }
    }
  }
}

/// Invokes `body(dd, ud, du, uu)` over the four corners of the 2x2 block the
/// two spins of `pair` span, for every setting of the remaining spins.
template <typename Body>
void forEachCouplingBlock(std::size_t nstates, SpinPair pair, Body body) {
  const std::size_t hiMask = spinMask(pair.high());
  const std::size_t loMask = spinMask(pair.low());
  const std::size_t hiBlock = hiMask << 1;
  const std::size_t loBlock = loMask << 1;
  for (std::size_t outer = 0; outer < nstates; outer += hiBlock) {
    for (std::size_t mid = outer; mid < outer + hiMask; mid += loBlock) {
      for (std::size_t dd = mid; dd < mid + loMask; ++dd) {
        body(dd, dd + hiMask, dd + loMask, dd + hiMask + loMask);
      }
    }
  }
}

}  // namespace

void axpySigmaX(double c, const StateVector& in, SpinIndex spin,
                StateVector& out) {
  checkOperands(in, out, spin);
  const double* i = in.data();
  double* o = out.data();
  forEachSpinPairIndex(in.size(), spin, [&](std::size_t lo, std::size_t hi) {
    const std::size_t l = 2 * lo;
    const std::size_t h = 2 * hi;
    o[h] += c * i[l];
    o[h + 1] += c * i[l + 1];
    o[l] += c * i[h];
    o[l + 1] += c * i[h + 1];
  });
}

void axpySigmaY(double c, const StateVector& in, SpinIndex spin,
                StateVector& out) {
  // sigma_y |down> = -i|up>, sigma_y |up> = +i|down>.
  checkOperands(in, out, spin);
  const double* i = in.data();
  double* o = out.data();
  forEachSpinPairIndex(in.size(), spin, [&](std::size_t lo, std::size_t hi) {
    const std::size_t l = 2 * lo;
    const std::size_t h = 2 * hi;
    o[h] += c * i[l + 1];
    o[h + 1] -= c * i[l];
    o[l] -= c * i[h + 1];
    o[l + 1] += c * i[h];
  });
}

void axpySigmaZ(double c, const StateVector& in, SpinIndex spin,
                StateVector& out) {
  // sigma_z = +1 on |up>, -1 on |down>.
  checkOperands(in, out, spin);
  const double* i = in.data();
  double* o = out.data();
  forEachSpinPairIndex(in.size(), spin, [&](std::size_t lo, std::size_t hi) {
    const std::size_t l = 2 * lo;
    const std::size_t h = 2 * hi;
    o[h] += c * i[h];
    o[h + 1] += c * i[h + 1];
    o[l] -= c * i[l];
    o[l + 1] -= c * i[l + 1];
  });
}

void applySigmaX(const StateVector& in, SpinIndex spin, StateVector& out) {
  checkOperands(in, out, spin);
  const double* i = in.data();
  double* o = out.data();
  forEachSpinPairIndex(in.size(), spin, [&](std::size_t lo, std::size_t hi) {
    const std::size_t l = 2 * lo;
    const std::size_t h = 2 * hi;
    o[h] = i[l];
    o[h + 1] = i[l + 1];
    o[l] = i[h];
    o[l + 1] = i[h + 1];
  });
}

void applySigmaY(const StateVector& in, SpinIndex spin, StateVector& out) {
  checkOperands(in, out, spin);
  const double* i = in.data();
  double* o = out.data();
  forEachSpinPairIndex(in.size(), spin, [&](std::size_t lo, std::size_t hi) {
    const std::size_t l = 2 * lo;
    const std::size_t h = 2 * hi;
    o[h] = i[l + 1];
    o[h + 1] = -i[l];
    o[l] = -i[h + 1];
    o[l + 1] = i[h];
  });
}

void applySigmaZ(const StateVector& in, SpinIndex spin, StateVector& out) {
  checkOperands(in, out, spin);
  const double* i = in.data();
  double* o = out.data();
  forEachSpinPairIndex(in.size(), spin, [&](std::size_t lo, std::size_t hi) {
    const std::size_t l = 2 * lo;
    const std::size_t h = 2 * hi;
    o[h] = i[h];
    o[h + 1] = i[h + 1];
    o[l] = -i[l];
    o[l + 1] = -i[l + 1];
  });
}

void axpyCoupling(double jx, double jy, double jz, const StateVector& in,
                  SpinPair pair, StateVector& out) {
  checkPair(in, out, pair);

  // sigma_x.sigma_x and sigma_y.sigma_y both connect |dd><->|uu| and
  // |ud><->|du|, with relative sign +1 and -1 respectively, so only these two
  // combinations ever appear.
  const double flipBoth = jx - jy;  // couples |dd> <-> |uu>
  const double flipOne = jx + jy;   // couples |ud> <-> |du>

  const double* i = in.data();
  double* o = out.data();
  forEachCouplingBlock(
      in.size(), pair,
      [&](std::size_t dd, std::size_t ud, std::size_t du, std::size_t uu) {
        const std::size_t d = 2 * dd;
        const std::size_t a = 2 * ud;
        const std::size_t b = 2 * du;
        const std::size_t u = 2 * uu;

        o[u] += flipBoth * i[d] + jz * i[u];
        o[u + 1] += flipBoth * i[d + 1] + jz * i[u + 1];
        o[d] += flipBoth * i[u] + jz * i[d];
        o[d + 1] += flipBoth * i[u + 1] + jz * i[d + 1];

        o[a] += flipOne * i[b] - jz * i[a];
        o[a + 1] += flipOne * i[b + 1] - jz * i[a + 1];
        o[b] += flipOne * i[a] - jz * i[b];
        o[b + 1] += flipOne * i[a + 1] - jz * i[b + 1];
      });
}

void axpyDiagonal(double c, std::span<const double> diagonal,
                  const StateVector& in, StateVector& out) {
  assert(diagonal.size() == in.size() && "diagonal does not match the space");
  assert(in.size() == out.size() && "state vectors span different spaces");
  const double* d = diagonal.data();
  const double* i = in.data();
  double* o = out.data();
  const std::size_t n = in.size();
  // One unit-stride pass, trivially vectorised -- this replaces every sigma_z
  // contribution in the Hamiltonian.
  for (std::size_t k = 0; k < n; ++k) {
    const double f = c * d[k];
    o[2 * k] += f * i[2 * k];
    o[2 * k + 1] += f * i[2 * k + 1];
  }
}

void axpyTransverse(double cx, double cy, const StateVector& in,
                    SpinIndex spin, StateVector& out) {
  checkOperands(in, out, spin);
  // In the (down, up) basis the combined operator is
  //     cx sigma_x + cy sigma_y = [[0, cx + i cy], [cx - i cy, 0]].
  //
  // This is the one kernel the interleaved layout does not suit: it mixes the
  // real and imaginary parts of each amplitude, so the compiler has to cross
  // lanes. Measured 3% slower than under the split layout. See StateVector for
  // why the layout was changed anyway, and why it did not pay off.
  const double* i = in.data();
  double* o = out.data();
  forEachSpinPairIndex(in.size(), spin, [&](std::size_t lo, std::size_t hi) {
    const std::size_t l = 2 * lo;
    const std::size_t h = 2 * hi;
    const double ar = i[l];
    const double ai = i[l + 1];
    const double br = i[h];
    const double bi = i[h + 1];
    o[l] += cx * br - cy * bi;
    o[l + 1] += cx * bi + cy * br;
    o[h] += cx * ar + cy * ai;
    o[h + 1] += cx * ai - cy * ar;
  });
}

void axpyFlipFlop(double c, const StateVector& in, SpinPair pair,
                  StateVector& out) {
  checkPair(in, out, pair);
  // Only the two anti-aligned corners are touched, so this moves half the data
  // the general coupling kernel does -- that part is a genuine win, and it is
  // where most of the speed in this file comes from.
  //
  // Two further attempts on this kernel measured well in isolation and did
  // nothing end to end: special-casing low == 0, whose innermost run is a
  // single element (1.86x on that pair alone), and interleaving the amplitudes
  // (1.7x on the kernel). Both are described in StateVector. The propagator is
  // limited by bytes moved, not by how well this loop issues.
  const double* i = in.data();
  double* o = out.data();
  forEachCouplingChannel(in.size(), pair, spinMask(pair.high()),
                         spinMask(pair.low()),
                         [&](std::size_t a, std::size_t b) {
                           o[a] += c * i[b];
                           o[a + 1] += c * i[b + 1];
                           o[b] += c * i[a];
                           o[b + 1] += c * i[a + 1];
                         });
}

void axpyFlipBoth(double c, const StateVector& in, SpinPair pair,
                  StateVector& out) {
  checkPair(in, out, pair);
  const double* i = in.data();
  double* o = out.data();
  forEachCouplingChannel(in.size(), pair, 0,
                         spinMask(pair.high()) + spinMask(pair.low()),
                         [&](std::size_t d, std::size_t u) {
                           o[u] += c * i[d];
                           o[u + 1] += c * i[d + 1];
                           o[d] += c * i[u];
                           o[d + 1] += c * i[u + 1];
                         });
}

}  // namespace spinsim
