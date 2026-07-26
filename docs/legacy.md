# The Fortran original

The `legacy/` directory holds the program this one replaces, kept verbatim. You
need this document for two reasons: the archived results in `../CPMGZ*` and
`../Rabi_*` were produced by it, and several of its defects change how those
results should be read.

## What it was

A Chebyshev-based spin-dynamics program, ~4500 lines of Fortran 77 across
thirteen files, written around 2019 and used for 2D-NMR echo modelling. The
physics was sound and is preserved here.

```
chebNMR2D-v6.f      main program, MPI master-worker over realisations
├── chparsgenPDDG   reads the @-script (@file, @CYCLE n … @ENDCYCLE)
│   └── chparsFPDDG parses one .dat file by substring matching
├── gendipnmr2D     random 2D geometry, Jz ∝ 1/r³
└── chebsdPDDG      runs the sequence for one realisation
    ├── chinitPDDG      initial state
    ├── chstepsPDDGnmr  real time, Chebyshev + rjbesl   ← the hot path
    ├── chimstepsPDDG   imaginary time, + ribesl
    └── chspecsteps     instantaneous gates             ← dead code, see below
```

## The problem it had

Not the physics — the fact that **the experiment was chosen by editing source and
copying the directory.** The tree contains 629 `.f` files across 47 directories,
of which only 29 are distinct. Four things varied, and each should have been
configuration:

| Varied | Where it was hard-coded |
|---|---|
| Which correlator to measure | `chebsdPDDGnmr.f:110` — `addHSx` in CPMGZ, `addHSz` in Rabi_10_SzSz |
| Which diagnostic files to write | `chebNMR2D-v6.f:342-377`, blocks commented in or out |
| Number of realisations | `chebNMR2D-v6.f:17` — `parameter (Nit=180)` |
| Spin count, τ, sequence | `.dat` files plus `maxLtot` in `chsdpar.h` |

Roughly forty of the forty-seven directories differ in nothing but the first.
The `Observable` registry is what collapsed them.

## Defects found, and what they mean for the archived data

These were found while porting. Each is verified.

**1. Error bars are the wrong quantity.** `chebNMR2D-v6.f:333` and `:425` print
`sqrt(M2/n)` under the heading of an uncertainty. That is the population standard
deviation — the spread between realisations — not the uncertainty of the mean,
which is smaller by √n. Any claim of statistical significance drawn from the old
output should be re-derived. This code reports both, separately named.

**2. Results depended on the MPI schedule.** `chebNMR2D-v6.f:468` does
`iRandT(ifile) = iRandT(ifile) + it*11489` on a worker-local copy that is never
reset. A realisation's random numbers therefore depended on how many jobs that
rank had already handled — that is, on scheduling. Two runs of the same input
could not be expected to agree. Compounding it, `srand` is never called anywhere,
so `rand()` always began from the default seed, and "seeding" meant discarding
values in a loop (`do i=1,lociRand; xx=rand()`), costing up to ~2×10⁶ calls per
realisation. Here a realisation is a pure function of (seed, index).

**3. A normalisation error in the π/2 gates.** `chspecsteps.f:277` and `:286` set
`ar = dsqrt(2.0D0)` where `exp(-i(π/4)σx)` requires 1/√2, so each application
doubled the norm. Harmless in practice — that file is dead code, see below — but
fixed here, and the gate unitarity test would catch it now.

**4. Memory sized by the compile-time maximum.** Every routine declared its
buffers as `psiR(0:maxStat-1)` with `maxStat = 2^24`, so a ten-spin run reserved
about 0.8 GB of locals to use 16 KB. Here the state is sized by the run.

**5. Nested loops rejected.** `chparsgenPDDG.f` returns status `-987` for a
nested `@CYCLE`, because the cycle was implemented as a pair of jump targets
rather than as structure. XY-8, XY-16 and concatenated DD were therefore never
run on this code. They are expressible here.

## Two things that look like bugs and are not

**The imperfect π pulse is deliberate.** `nmrtst2eP1.dat` sets Hx = 10700 with
τ = 0.000314159265, so Hx·τ = 3.3615 rad = **1.07π** — a 7% over-rotation.
That is the pulse-error physics the work was about, matching the Lang et al.
paper in `../Papers/`. Do not "fix" it.

**Time runs as exp(+iHt).** `chstepsPDDGnmr.f:175` sets the scaled Hamiltonian
to `−H/emax`, which makes the expansion realise exp(+iHt) rather than the
textbook exp(−iHt). Verified against exact diagonalisation and preserved here so
that comparison with the archived results stays valid. Every observable in use is
invariant under it. See [physics.md](physics.md#conventions).

## Dead code

`chspecstepsPDDGnmr.f` — 337 lines implementing fifteen instantaneous gates
selected by magic number — is **never reached**. Its only entry point is the
`SPECIAL PROPAGATION` keyword, which appears in none of the `.dat` files.

The real pulses are finite-duration: evolution under a large transverse field.
That is the physically interesting case, because the couplings keep acting during
the pulse. This code offers both, so the two can be compared by subtraction.

## Migrating input files

```bash
python3 tools/dat2toml.py ../CPMGZ/CPMG8Z/nmrCPMG8.dat > configs/cpmg8.toml
./build/spinsim validate configs/cpmg8.toml
```

All 134 archived `@`-scripts convert and validate. Three settings cannot be
recovered, because they lived in Fortran source rather than in any data file, and
the converter marks each with a `TODO`:

- the measured correlator (`chebsdPDDGnmr.f:110`),
- the realisation count (`chebNMR2D-v6.f:17`),
- the disorder strength `HzScale` (`chebNMR2D-v6.f:56`).

The converter is Python because it is a one-time tool that must be forgiving of
input the real parser should reject.

### Why the old format needed replacing

`chparsFPDDG` matched keywords by substring: `index(curline,'END')` closed a
block on any line containing those three letters anywhere. Errors were magic
integers — `-11`, `-83`, `-1701`, `-4444` — written to a log file the user had to
find, after which parsing continued. The TOML parser validates against a schema
and names the offending key.

## Building the original

Only needed to regenerate Bessel reference data. The Makefile points at a
cluster path (`/opt/ud/openmpi-1.8.8/bin/mpifort`) that no longer exists;
`gfortran` works for the standalone numerical routines:

```bash
cd tests/data
gfortran -O2 -std=legacy gen_bessel_reference.f ../../legacy/{rjbesl,ribesl}.f -o gen_bessel_reference
```

Building the full program single-rank would need an MPI stub. It has not been
done, because verification against exact diagonalisation is a stronger check than
agreement with the original would be.
