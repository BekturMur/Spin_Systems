# Legacy archive

This directory preserves the Fortran implementation, historical experiment
inputs, and research notes from which the modern C++ `spinsim` was developed.
It is archival material: the supported simulator, build system, tests, and
documentation are in the repository root.

## Contents

| Path | Contents |
|---|---|
| `experiments/cpmg/` | 12 historical CPMG experiment directories and their lightweight `.dat` inputs |
| `experiments/rabi/` | 33 Rabi and correlator experiment directories |
| `reference/` | Canonical source set used as the stable Fortran reference |
| `variants/` | Distinct historical source variants, named by the first 12 characters of their SHA-256 digest |
| `source-manifest.csv` | Mapping from all 779 pre-cleanup source paths to 45 canonical files |
| `materialize.py` | Verification and experiment reconstruction tool |
| `prototypes/` | Early Chebyshev implementations and root-level prototypes |
| `notes/` | LyX, TeX, BibTeX, and Mathematica research notes |
| `literature/` | Bibliographic note for third-party papers, which are not redistributed |

Exact duplicate Fortran, header, makefile, and plotting-script copies were
removed. Experiment data files were deliberately left in place: their directory
context and relative includes are part of the historical experiment definition.

Generated executables, simulation output, TeX build products, and local copies
of papers are not versioned. The omitted local results occupy several gigabytes
and can contain files above GitHub's size limit.

## Reconstruct an experiment

Run the tool from the repository root and choose one directory directly below
`experiments/cpmg/` or `experiments/rabi/`:

```bash
python3 legacy/materialize.py legacy/experiments/rabi/Rabi_10_SzSz \
  --out /tmp/Rabi_10_SzSz
```

The output directory must not already exist. Only tracked experiment artifacts
are copied, so unversioned multi-gigabyte result files are not pulled into the
reconstruction. When using a GitHub source archive without `.git`, the clean
files present in that archive are copied instead.

Verify every canonical path and checksum with:

```bash
python3 legacy/materialize.py --verify
```

The manifest keeps paths as they appeared before the repository cleanup. This
makes old notes and Git history searchable even though the physical duplicates
are gone.

## Scientific caveats

The archive is evidence and migration material, not the recommended simulator.
In particular, the original program had schedule-dependent random streams and
reported population spread as though it were uncertainty of the mean. It also
uses the convention `exp(+iHt)`. Read [`docs/legacy.md`](../docs/legacy.md) and
[`docs/physics.md`](../docs/physics.md) before comparing archived numbers.
