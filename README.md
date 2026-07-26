# spinsim

[![CI](https://github.com/BekturMur/Spin_Systems/actions/workflows/ci.yml/badge.svg)](https://github.com/BekturMur/Spin_Systems/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Exact quantum dynamics of interacting spin-1/2 ensembles, for modelling NMR echo
experiments and dynamical-decoupling sequences.

Given a set of spins, their couplings, and a pulse sequence, it propagates the
Schrödinger equation exactly — no master equation, no semi-classical bath, no
truncation of the Hilbert space — and reports echo signals averaged over
disorder. The practical ceiling is 16 spins comfortably, 18 with patience.

It is a rewrite of a Fortran 77 program used for this physics around 2019. See
[docs/legacy.md](docs/legacy.md) for what changed and why, including several
bugs found in the original that affect how its archived results should be read.

## Five minutes to a first result

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build            # 86 tests, should all pass

./build/spinsim describe configs/cpmg8z.toml
./build/spinsim run configs/cpmg8z.toml --realizations 20 --out /tmp/echo.csv
```

That last command runs a CPMG echo experiment on ten dipolar-coupled spins and
writes a CSV whose header names every column. It takes about a fifth of a second.

Needs CMake ≥ 3.28 and a C++23 compiler. On macOS, `brew install cmake`; LAPACK
comes from the system Accelerate framework. Catch2 and toml++ are fetched
automatically at configure time, so the first build needs network access.
`gfortran` is only needed to regenerate the Bessel reference data.

On Debian or Ubuntu, install the system dependencies with
`sudo apt install cmake ninja-build g++ liblapack-dev libblas-dev`.

## Where to read next

| Document | What it covers |
|---|---|
| [docs/physics.md](docs/physics.md) | The Hamiltonian, the Chebyshev method, what the observables mean, and **the sign and ordering conventions** — read this before trusting any number |
| [docs/configuration.md](docs/configuration.md) | Complete TOML reference with worked examples, including XY-8 |
| [docs/architecture.md](docs/architecture.md) | How the code is layered, and how to add an observable, a step kind, or a geometry |
| [docs/legacy.md](docs/legacy.md) | The Fortran original: what it did, what was wrong with it, how to migrate old input files |
| [docs/claude-science-prompt.md](docs/claude-science-prompt.md) | A briefing for planning research with this tool |

If you are here to run experiments rather than change code, physics and
configuration are enough.

## What it can and cannot do

Can:

- Arbitrary anisotropic Heisenberg couplings and per-spin local fields.
- Pulse sequences as a nested program — `evolve`, `pulse`, `measure`, `repeat` —
  so composite sequences like XY-8 and XY-16 are expressible.
- Both realistic finite-duration pulses (evolution under a strong field) and
  ideal instantaneous rotations, which lets you separate error sources by
  subtraction.
- Monte-Carlo averaging over randomly generated 2D dipolar geometries and
  Gaussian field disorder, with correct error bars.
- Real-time and imaginary-time propagation.

Cannot, today:

- Exceed about 18 spins. The state is 2^L amplitudes; this is the binding limit.
- Report entanglement entropy or any subsystem quantity — there is no reduced
  density matrix.
- Take an explicit coupling table from a config file, so chains, lattices and
  custom topologies are out of reach without a small extension. The
  `Hamiltonian` class supports them; only the parser does not.
- Model open-system dynamics. Everything is unitary; decoherence has to come
  from the spin bath itself.

[docs/architecture.md](docs/architecture.md) estimates what each of those costs
to add. The first two are small.

## Why the numbers can be trusted

The propagator is checked against an entirely different algorithm — dense
eigendecomposition via LAPACK — and agrees to better than 1e-10. Norm and energy
are conserved to 1e-11 over 10⁴ steps with autonormalisation switched off, so
expansion error is measured rather than hidden. The Bessel coefficients are
checked against the original Cody routines compiled from Fortran, and the sparse
kernels against dense Kronecker-product Pauli matrices built independently.

Results are also bit-for-bit reproducible: a run on one thread and a run on
eighteen produce byte-identical files, because each Monte-Carlo realisation
draws from a stream fixed by (seed, index) and partial results are combined in
index order. The original could not do this.

Details in [docs/physics.md](docs/physics.md#verification).

## Layout

```
include/spinsim/, src/    library, split into core / sequence / ensemble / io
apps/spinsim/             command-line front end
tests/                    86 tests, plus bench_propagator
tools/dat2toml.py         converter for the legacy .dat input format
configs/                  worked examples
docs/                     physics, configuration, architecture, migration
legacy/                   categorized Fortran archive and experiment inputs
```

About 3500 lines of library and application code, 2100 lines of tests.

## Legacy archive

The former experiment-specific Fortran trees are under [`legacy/`](legacy/).
Their 779 historical source paths contained only 45 distinct files, so exact
copies now live once in `legacy/reference/` or `legacy/variants/`. The
content-addressed [`source-manifest.csv`](legacy/source-manifest.csv) records
the original path, canonical path, and SHA-256 checksum for every copy.

To reconstruct an old experiment as a self-contained directory:

```bash
python3 legacy/materialize.py legacy/experiments/rabi/Rabi_10_SzSz \
  --out /tmp/Rabi_10_SzSz
python3 legacy/materialize.py --verify
```

Simulation outputs and third-party papers are intentionally excluded from Git.
The tracked inputs, notes, canonical sources, and reconstruction metadata remain.
See the [archive guide](legacy/README.md) before interpreting old results.

## Contributing and license

Contributions are welcome; start with [CONTRIBUTING.md](CONTRIBUTING.md). This
project is released under the [MIT License](LICENSE). Historical third-party
material retained for provenance is identified in [THIRD_PARTY.md](THIRD_PARTY.md).
