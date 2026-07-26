#pragma once

#include <span>

#include "spinsim/core/Bits.hpp"
#include "spinsim/core/StateVector.hpp"

/// Pauli-operator kernels: ports of addHSx/addHSy/addHSz, addJSxyz and
/// actSx/actSy/actSz from chstepsPDDGnmr.f.
///
/// Convention: these apply the **Pauli** matrices sigma, not the spin operators
/// S = sigma/2. The legacy code did the same and folded the factor of two into
/// the coefficients it passed (see chstepsPDDGnmr.f:175, "= -0.5*H*2"). Keeping
/// the kernels in sigma and the physics in S means exactly one place --
/// Hamiltonian::applyTwiceNormalized -- has to know about the factor.
///
/// Every function reads `in` and accumulates into (axpy*) or overwrites
/// (apply*) `out`. `in` and `out` must not alias; the legacy code relied on
/// this silently, here it is asserted.
namespace spinsim {

/// out += c * sigma_x(spin) * in
void axpySigmaX(double c, const StateVector& in, SpinIndex spin,
                StateVector& out);

/// out += c * sigma_y(spin) * in
void axpySigmaY(double c, const StateVector& in, SpinIndex spin,
                StateVector& out);

/// out += c * sigma_z(spin) * in
void axpySigmaZ(double c, const StateVector& in, SpinIndex spin,
                StateVector& out);

/// out = sigma_x(spin) * in
void applySigmaX(const StateVector& in, SpinIndex spin, StateVector& out);

/// out = sigma_y(spin) * in
void applySigmaY(const StateVector& in, SpinIndex spin, StateVector& out);

/// out = sigma_z(spin) * in
void applySigmaZ(const StateVector& in, SpinIndex spin, StateVector& out);

/// out += (jx sx.sx + jy sy.sy + jz sz.sz) * in, for the two spins in `pair`.
///
/// Port of addJSxyz. The three products share one traversal of the state
/// vector, which is why the legacy code fused them (addJSxxSyy and addJSzz are
/// the unfused versions, left behind and unused).
void axpyCoupling(double jx, double jy, double jz, const StateVector& in,
                  SpinPair pair, StateVector& out);

// -- specialised kernels -------------------------------------------------
//
// The general kernels above are what the tests measure correctness against.
// These are the ones the propagator actually calls, and they exist because the
// generic form does far more work than the physics needs:
//
//  * every sigma_z term -- single-spin and coupling alike -- is diagonal, so
//    all of them collapse into one precomputed vector and a single pass;
//  * the secular dipolar coupling has jx == jy, which makes the
//    |down down> <-> |up up> channel vanish identically, halving what is left.
//
// The legacy code exploited neither: chstepsPDDGnmr.f called addHSx, addHSy and
// addHSz separately per spin and addJSxyz per pair, so a 14-spin run made 105
// full passes over the state vector for every order of the expansion.

/// out += c * diag[k] * in[k], elementwise.
void axpyDiagonal(double c, std::span<const double> diagonal,
                  const StateVector& in, StateVector& out);

/// out += (cx sigma_x + cy sigma_y)(spin) * in, in one traversal.
void axpyTransverse(double cx, double cy, const StateVector& in,
                    SpinIndex spin, StateVector& out);

/// out += c * (|up down> <-> |down up>) for `pair` -- the flip-flop channel,
/// which is (jx + jy)/2 times the exchange of the two anti-aligned states.
void axpyFlipFlop(double c, const StateVector& in, SpinPair pair,
                  StateVector& out);

/// out += c * (|down down> <-> |up up>) for `pair`. Zero whenever jx == jy, so
/// the propagator skips it for dipolar couplings.
void axpyFlipBoth(double c, const StateVector& in, SpinPair pair,
                  StateVector& out);

// Measured and rejected: fusing every coupling into a single traversal, with
// the output summed in registers and the partner amplitudes read by XOR, is
// twice as slow as the per-pair kernels above (1054 ms against 498 ms for
// twenty steps at fourteen spins). A pair contributes only when its two spins
// are anti-aligned, so half the gathered reads are wasted, and scattered reads
// cost more than the short sequential runs they were meant to replace. The
// per-pair loop nest stays.

}  // namespace spinsim
