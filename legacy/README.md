# Spin Systems

Research code and archived experiment definitions for exact quantum dynamics of
interacting spin-1/2 ensembles, with an emphasis on NMR echoes, CPMG sequences,
Rabi oscillations, and dynamical decoupling.

The actively maintained implementation is [`spinsim/`](spinsim/README.md), a
C++23 rewrite of the original Fortran 77 simulator. It replaces dozens of
experiment-specific source copies with one tested executable driven by TOML
configuration files.

## Repository map

| Path | Contents |
|---|---|
| [`spinsim/`](spinsim/) | C++23 library and CLI, 86 tests, example TOML configurations, migration tool, documentation, and a verbatim reference copy of the main legacy sources |
| `CPMGZ*` | Archived CPMG experiment variants: Fortran sources and lightweight `.dat` inputs |
| `Rabi*` | Archived Rabi/correlator variants for different spin counts, interactions, and pulse parameters |
| `Cheb/` and root `*.f` | Earlier Chebyshev-propagation sources and numerical routines |
| `Graphics_CPMG/` and `graphika.py` | Plotting and exploratory analysis scripts/notebooks |
| root `*.lyx`, `*.tex`, `*.bib`, `*.nb` | Research notes and derivations |

Generated executables, build trees, simulation output, TeX products, and local
copies of third-party papers are intentionally not versioned. The local output
archive is about 4.6 GiB and contains individual files above GitHub's 100 MiB
limit; all source code and lightweight experiment inputs needed to reproduce
new results are retained.

## Build and verify

From the repository root:

```bash
cmake -S spinsim -B spinsim/build -DCMAKE_BUILD_TYPE=Release
cmake --build spinsim/build -j
ctest --test-dir spinsim/build --output-on-failure
```

Requirements are CMake 3.28 or newer, a C++23 compiler, and LAPACK. The first
configure needs network access to fetch pinned versions of Catch2 and toml++.

Run a small CPMG experiment:

```bash
./spinsim/build/spinsim describe spinsim/configs/cpmg8z.toml
./spinsim/build/spinsim run spinsim/configs/cpmg8z.toml \
  --realizations 20 --out results/cpmg8z.csv
```

## Scientific and numerical notes

- The state vector contains `2^L` complex amplitudes; 16 spins is comfortable
  and 18 is practical with patience.
- The Chebyshev propagator is verified against dense LAPACK
  eigendecomposition, while Bessel coefficients are checked against the legacy
  Cody Fortran routines.
- Monte-Carlo realisations are reproducible across thread counts because every
  stream is derived from `(seed, realisation index)` and results are accumulated
  in a fixed order.
- For compatibility with the archived calculations, real-time propagation uses
  `exp(+iHt)`, not the textbook `exp(-iHt)`. See
  [`spinsim/docs/physics.md`](spinsim/docs/physics.md) before comparing signs or
  phases with an external calculation.
- The original archive has schedule-dependent random streams and mislabelled
  error bars. See [`spinsim/docs/legacy.md`](spinsim/docs/legacy.md) before using
  old `.out` files quantitatively.

## Documentation

- [Physics and conventions](spinsim/docs/physics.md)
- [Configuration reference](spinsim/docs/configuration.md)
- [Architecture](spinsim/docs/architecture.md)
- [Legacy implementation and migration](spinsim/docs/legacy.md)

No license has been selected yet; absence of a license means reuse permission is
not granted automatically.

