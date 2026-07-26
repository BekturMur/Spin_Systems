# Physics and conventions

Read the conventions section before trusting any number this code produces. Two
of them — the sign of time evolution and the meaning of the correlator — will
mislead you if you assume the textbook version.

## The model

L spin-1/2 particles, state vector of 2^L complex amplitudes, evolving under

    H = Σ_i (Hx_i Sx_i + Hy_i Sy_i + Hz_i Sz_i)
      + Σ_{i<j} (Jx_ij Sx_i Sx_j + Jy_ij Sy_i Sy_j + Jz_ij Sz_i Sz_j)

with **S = σ/2**, so a field H_x drives Rabi flopping at angular frequency H_x.
Fields and couplings carry whatever units you choose; only their ratios and
their products with time matter.

Fields are piecewise constant. They hold one value for the whole of a sequence
step and may change at the next. There is no continuous time dependence inside a
step, so a resonant oscillating drive has to be treated in a rotating frame, or
approximated by subdividing into many short steps.

## Conventions

**Basis ordering.** Spin *s* owns bit *s* of the basis index, counting from the
least significant. A set bit means "up". So for three spins, index 5 = binary
101 is |up, down, up⟩ with spin 0 up. Spin indices are 0-based in the code and
**1-based in configuration files**, matching the legacy `.dat` format;
conversion happens at the parser.

**Pauli, not spin, inside the kernels.** `SpinOperators` applies σ, not S = σ/2.
The factors of one half live in `Hamiltonian::accumulate`, which is the single
place that knows the physical normalisation. If you write a new kernel, follow
that: kernels in σ, physics in S.

**Time evolution runs as exp(+iHt), not exp(−iHt).** This is inherited from the
Fortran original, where `chstepsPDDGnmr.f:175` sets the scaled Hamiltonian to
`−H/emax`; the expansion then sums `(−i)ⁿ Tₙ(−H/emax) = (+i)ⁿ Tₙ(H/emax)`. It is
verified against exact diagonalisation and documented rather than silently
changed, because changing it would break comparison with the archived results.

Equivalently: the code treats the tabulated fields and couplings as −H. Every
observable currently implemented is invariant under this — they are built from
correlators that are even under time reversal — but it matters the moment you
compare against an analytic formula, or add an observable sensitive to the
direction of precession.

**Secular truncation is baked into the geometry generator.** The dipolar
generator emits Jx = Jy = −Jz/2, the secular (high-field) form. If you want a
different coupling anisotropy you must set couplings directly through the
`Hamiltonian` API; the config file cannot express it yet.

## Propagation

The step operator is expanded in Chebyshev polynomials of the rescaled
Hamiltonian:

    exp(-i α Ĥ) = J₀(α) + 2 Σ_{n≥1} (-i)ⁿ Jₙ(α) Tₙ(Ĥ),   Ĥ = -H/E,  α = E·τ

where E is an upper bound on ‖H‖ so that Ĥ has spectrum inside [−1, 1], and Jₙ
are Bessel functions of the first kind. Imaginary time uses the same structure
with modified Bessel functions Iₙ and alternating signs, giving exp(+Hτ).

Three things are worth knowing:

**The bound is the sum of term norms**, Σ_i |H_i|/2 + Σ_ij |J_ij|/4, not the true
spectral radius. Chebyshev needs a genuine upper bound — underestimating it
destroys convergence — and this one is cheap. It is also tight: measured against
exact diagonalisation it overestimates by only 7–11% for the disordered dipolar
systems used here, so there is nothing to gain from a better estimate.

**The order is chosen from the coefficients, not from a formula.** The expansion
keeps terms until three consecutive Bessel coefficients are all below `epsilon`
*and* decreasing. Three, because one or two can dip by accident near a zero of
Jₙ. The default `epsilon` is 1e-7.

**The order grows linearly with α = E·τ.** A step ten times longer costs roughly
ten times more. If a run is slow, look at the step length before anything else;
`bench_propagator` prints the order it selected.

If the coefficients have not decayed by the maximum order, the run stops with a
`ChebyshevConvergenceError` naming the bound, the step, and α. Shorten the step
or raise `epsilon`.

## Observables

**`corr:<axis>` — the echo signal, and the main quantity.**

At t = 0 a companion state is prepared as A|ψ(0)⟩ with A = Σ_i σ_b(i), where b is
`measurement.reference_axis`. Both states are then propagated through the *same*
sequence. At each measurement point the observable reports

    ⟨ref(t)| Σ_i σ_a(i) |ψ(t)⟩

as real and imaginary parts, where a is the axis in the observable name. On a
random typical initial state this estimates the infinite-temperature two-time
correlator Tr[A(0) B(t)] / Tr[1] — the NMR echo signal.

Two consequences worth internalising:

- At t = 0 with a = b, the value is ‖A|ψ⟩‖² ≈ L. Seeing roughly the spin count in
  the first row is the sanity check that the run is set up correctly.
- The reference axis and the measured axis are independent. `reference_axis = "x"`
  with `observables = ["corr:x"]` reproduces the old CPMGZ runs; `"z"` with
  `corr:z` reproduces Rabi_10_SzSz. In the Fortran these were separate compiled
  binaries in separate directories.

**`mag:<spin>:<axis>`** — ⟨ψ(t)| σ_a(s) |ψ(t)⟩ for one spin, 1-based in config.

**`norm`** — norms of both states. Under real-time evolution these are conserved
exactly, so any drift is accumulated expansion error, reported directly. Include
it in production runs; it costs nothing and it is the cheapest possible check
that a result is numerically sound.

**`energy`** — ⟨ψ|H|ψ⟩, also conserved, and an independent check on the same
thing. The original recorded neither.

## Typicality

The default initial state is a random vector drawn from the unit sphere in the
full 2^L space. Such a state reproduces infinite-temperature trace averages to
relative accuracy O(2^(−L/2)) — about 0.4% at L = 16 — which is why a single
realisation is already meaningful and why the ensemble average converges quickly.

This is a physics assumption, not a numerical one. If you need a finite
temperature or a specific prepared state, typicality does not apply and the
initial state must be constructed deliberately.

## Disorder and geometry

Each Monte-Carlo realisation regenerates, in this order:

1. **Positions.** Spin 0 at the origin; the rest uniform in a square of
   half-width √L/2. Coupling Jz = `dipolar_scale`/r³, with Jx = Jy = −Jz/2.
2. **Field disorder.** A Gaussian offset of width `field_disorder` added to each
   Hz_i.
3. **The initial state.**

The order matters: it fixes the random stream, and changing it changes every
result even at the same seed.

**A caveat inherited from the original.** A candidate position is rejected only
if it falls within 0.05 of an earlier spin in *both* coordinates at once. That is
not a minimum-distance test — two spins can be arbitrarily close along a
diagonal — so occasional very large couplings occur. Measured at L = 16 over 40
realisations, the largest |Jz| was 4.4×10³ against a typical scale of order 1.
Such a realisation inflates the spectral bound and therefore the expansion order,
which is why run times fluctuate between realisations. It was kept as-is to stay
faithful; if it distorts a result, that is a physics decision to revisit.

## Verification

What each test actually establishes:

| Check | Test tag | Establishes |
|---|---|---|
| Against `zheev` matrix exponential, L ≤ 6 | `[propagator][exact]` | The whole propagation chain is correct to 1e-10 |
| Norm and energy over 10⁴ steps, autonorm off | `[propagator][conservation]` | Expansion error is bounded and measured, not masked |
| Bessel against `rjbesl.f`/`ribesl.f` | `[bessel]` | 400 orders at 14 arguments agree to 1e-14 |
| Kernels against dense Pauli matrices, L ≤ 5 | `[kernels]` | Bit arithmetic matches textbook operators |
| σx σy = i σz on the kernels | `[algebra]` | Fixes the relative sign of σy without the dense reference |
| Single spin in a transverse field | `[propagator]` | Reproduces analytic Rabi precession |
| 1 thread versus 18 | `[determinism]` | Byte-identical output |

Run a subset by tag:

```bash
./build/spinsim_tests "[propagator][exact]"
```

The dense reference in `tests/DenseReference.cpp` is built from explicit
Kronecker products and never reuses the bit arithmetic under test, so agreement
between the two is evidence rather than tautology. Keep it that way.
