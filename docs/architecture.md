# Architecture

## The organising rule

`core/` knows nothing about files, threads, or configuration. It is pure
numerics: given a Hamiltonian and a state, evolve it. Everything about *how a
study is organised* — parsing, averaging, parallelism, output — sits in the
layers above and depends downward only.

This is the fix for what went wrong in the original. There, MPI calls, file I/O,
and physics were interleaved in the same routines, so the 2D-NMR modification had
to comment out half the broadcasts once the couplings became worker-local. Here,
adding a cluster backend touches one file and no physics.

```
apps/spinsim          command line
   │
io/                   TOML parsing, experiment assembly, CSV writing
   │
ensemble/             RNG, geometry, disorder, averaging, thread pool
   │
sequence/             the pulse-sequence program and its interpreter
   │
core/                 state, operators, Hamiltonian, Chebyshev, gates, observables
```

## core/

| File | Responsibility |
|---|---|
| `Bits.hpp` | Spin indices and masks. `SpinPair` carries the invariant `low < high` in the type — the original stated it only in a comment and relied on every call site. |
| `StateVector` | 2^L complex amplitudes, interleaved (`data()[2k]`, `data()[2k+1]`). Sized by the run, not by a compile-time maximum. |
| `SpinOperators` | The hot kernels. Pauli matrices, not spin operators. |
| `Hamiltonian` | Fields, couplings, spectral bound, and application to a state. Owns the cached diagonal. |
| `Bessel` | Jₙ and Iₙ sequences by Miller's backward recurrence. |
| `Chebyshev` | Expansion order selection and the coefficient table. |
| `Propagator` | The recurrence. Owns its scratch buffers. |
| `Gate` | Instantaneous rotations and two-spin unitaries. |
| `Observable` | The measurement registry. |

### Why the kernels look the way they do

`Hamiltonian::accumulate` dispatches to three kinds of kernel, and the split is
physics, not micro-optimisation:

- **Everything diagonal in one pass.** Every σ_z term — local fields and the
  Ising part of every coupling — is diagonal, so all of them collapse into one
  cached vector and a single streaming pass. Rebuilding it costs about one
  Hamiltonian application and is reused across all ~60 expansion orders in a
  step. Invalidated by any change to a field or coupling.
- **Transverse fields fused per spin.** σ_x and σ_y for one spin share a
  traversal.
- **Coupling split into two channels.** σxσx and σyσy connect |↑↓⟩↔|↓↑⟩ with
  weight (Jx+Jy) and |↓↓⟩↔|↑↑⟩ with (Jx−Jy). For the secular dipolar form
  Jx = Jy the second channel **vanishes identically**, so it is skipped: four
  operations per block instead of sixteen, touching two corners instead of four.

Together these are 3.4× faster than the direct transcription of the Fortran,
which called `addHSx`, `addHSy`, `addHSz` separately per spin and `addJSxyz` per
pair — 105 full passes over the state vector for every expansion order at 14
spins.

### What has been tried and rejected

The kernels are memory-bandwidth-bound above L ≈ 15, and several plausible ideas
do not help. Each is documented at the point in the code where someone would be
tempted to try it again. Do not repeat them without new evidence:

| Idea | Result |
|---|---|
| Fuse all couplings into one traversal, summing in registers | **2× slower.** Half the gathered reads are wasted (a pair contributes only when its spins are anti-aligned) and scattered reads cost more than short sequential runs. |
| Process both propagated states in one kernel pass | **0.6×.** Doubles the working set; it falls out of cache. |
| Special-case the pair involving spin 0, whose inner run is one element | 1.86× on that pair alone, **no end-to-end change.** In a real all-to-all system those pairs come with small partners and small blocks. |
| Tighter spectral bound to reduce expansion order | **Closed off.** The sum-of-norms bound overestimates by only 7–11%. |
| Interleaved instead of split real/imaginary storage | 1.7× on the dominant kernel in isolation, **no end-to-end change.** Kept because it is equivalent and simpler, not because it is faster. |
| `-mcpu=native`, FMA contraction | ~1%. Not worth losing reproducibility across optimisation levels. |
| GPU | **Not possible here.** Apple Silicon has unified memory, so the GPU shares the same bandwidth; and Metal has no FP64, which would destroy the verification chain. |

The lesson, learned the hard way more than once: **an isolated kernel benchmark
does not predict end-to-end behaviour.** The propagator keeps three state vectors
and the diagonal live, close to a megabyte at L = 16, and is limited by bytes
moved. Measure `bench_propagator`, not a loop.

Machine ceiling, measured by STREAM triad: 119.7 GB/s on one thread, 237.7 GB/s
on eighteen. The simulation reaches about 58 GB/s single-threaded at L = 16, so
threads scale sublinearly and the memory system is the limit.

## sequence/

The pulse sequence is an AST:

```cpp
using Step = std::variant<Evolve, Pulse, Measure, Repeat>;
struct Repeat { std::size_t count; std::vector<Step> body; };
```

`Repeat` nests, which is the whole point — the original implemented `@CYCLE` as a
pair of jump targets (`imark`) and returned error `-987` for a nested one, so
XY-8 and concatenated sequences were unreachable.

`Sequence::run` walks the tree, propagating **both** the main state and the
reference companion through every `Evolve` and `Pulse`, which is what makes the
two-time correlator meaningful. Field overrides are installed for the duration of
a step and restored after.

## ensemble/

| File | Responsibility |
|---|---|
| `Rng` | One stream per realisation, derived from (seed, index) by SplitMix64. |
| `Ensemble` | Initial and reference states, dipolar geometry, field disorder, and the `Accumulator`. |
| `EnsembleRunner` | The interface, plus `LocalRunner` — a `std::jthread` pool. |

Two properties matter more than they look:

**Reproducibility.** A realisation is a pure function of (config, index). The
original derived its randomness from a counter that accumulated in a worker-local
variable and was never reset, so a realisation's numbers depended on how many
jobs that MPI rank had already run — results depended on scheduling.

**Order-independent averaging.** Realisations run in whatever order threads pick
them up, but their outputs are folded into the statistics strictly in index
order. Floating-point summation is not associative, so without that discipline
the answer would drift with the thread count. `Accumulator` also implements
Chan's parallel merge, used by the tests and available for a future distributed
backend.

## io/

`Config` parses TOML into a `RunConfig` and validates it, throwing `ConfigError`
with the offending key. `Experiment` assembles a realisation and runs it;
`writeCsv` emits columns named by the observables themselves.

The original's output was a flat `double` array with a parallel array of record
lengths, and its layout existed only in the reader's head.

## Extension points

### Add an observable

1. Subclass `Observable` in `core/Observable.hpp`, implementing `name()`,
   `columns()` and `measure()`. Read what you need from `MeasureContext` —
   both states, the Hamiltonian, and the subsystem size are there.
2. Register the spelling in `makeObservable()`.
3. Test it against the dense reference in `tests/DenseReference.cpp`.

Nothing else changes: the CSV header comes from `columns()`, so output is
self-describing automatically. This is the abstraction that collapsed 47
near-identical Fortran directories into one binary.

### Add a step kind

Add a struct, add it to the `Step` variant, handle it in `runStep`, parse it in
`Config.cpp`. The variant makes the compiler find every place that needs
updating.

### Extensions the code is missing, with estimates

| Want | Cost | Where |
|---|---|---|
| Explicit J_ij in config — unlocks chains, lattices, XXZ | ~40 lines | `Config.cpp`; `Hamiltonian` already supports it |
| Reduced density matrix, entanglement entropy | ~80 lines | New `Observable`; `MeasureContext::subsystemSize` exists but is unused |
| 3D geometry | ~30 lines | `generateDipolarLattice2D`, alongside it |
| Parameter sweep driver | ~100 lines | New subcommand in `apps/` |
| Distributed backend | moderate | New `EnsembleRunner`; `core/` untouched |
| Single precision | large | Would roughly double throughput by halving bytes moved, but must be *shown* adequate against double, not assumed |

## Tests

86 tests. `tests/DenseReference.cpp` builds Pauli operators and matrix
exponentials from explicit Kronecker products and LAPACK, never reusing the bit
arithmetic under test — keep that independence or the tests stop being evidence.

```bash
ctest --test-dir build --output-on-failure
./build/spinsim_tests "[propagator][exact]"     # by tag
./build/bench_propagator 14 20                  # timings, not a test
```

Regenerating the Bessel reference data needs gfortran and is only necessary if
`legacy/reference/rjbesl.f` changes:

```bash
cd tests/data
gfortran -O2 -std=legacy gen_bessel_reference.f \
  ../../legacy/reference/{rjbesl,ribesl}.f -o gen_bessel_reference
./gen_bessel_reference > bessel_reference.txt
```

## Build

C++23, CMake ≥ 3.28. Catch2 and toml++ arrive via `FetchContent`, LAPACK from
Accelerate on macOS or `find_package` elsewhere.

`-ffp-contract=off` is deliberate: it keeps fused multiply-add out of the
kernels, so results do not shift between optimisation levels. It was measured to
cost about 1%, which is worth paying for reproducibility.
