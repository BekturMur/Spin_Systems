# Configuration reference

A run is one TOML file. Nothing else — no recompilation, no environment
variables. Check a file before running it:

```bash
./build/spinsim validate configs/my_run.toml
./build/spinsim describe configs/my_run.toml
```

Errors name the offending key and the file:

```
spinsim: configs/my_run.toml: 'measurement.observables' unknown observable
'corr:w'; expected one of corr:<axis>, mag:<spin>:<axis>, norm, energy
```

## A minimal file

Everything not listed has a default; these four sections are the minimum that
does something.

```toml
[system]
spins = 3

[measurement]
observables = ["corr:z"]

[[step]]
kind = "measure"

[[step]]
kind = "evolve"
tau = 0.1

[[step]]
kind = "measure"
```

## `[system]`

| Key | Default | Meaning |
|---|---|---|
| `spins` | required | Number of spins L. State is 2^L amplitudes; 16 is comfortable, 18 slow, 20 overnight. |
| `initial_state` | `"random"` | `random` (typical state — the physically relevant choice), `up`, `down`, or `basis`. |
| `basis_index` | `0` | Which basis state, when `initial_state = "basis"`. |

## `[field]`

Baseline local fields, in effect unless a step overrides them. Each component is
either one number applied to every spin, or an array with one entry per spin.

```toml
[field]
hx = 0.0
hz = [1.0, 1.0, 1.5]     # per-spin, length must equal `spins`
```

Omitted components are zero.

## `[ensemble]`

| Key | Default | Meaning |
|---|---|---|
| `realizations` | `1` | Monte-Carlo samples to average. |
| `seed` | `0` | Base seed. Same seed gives byte-identical output regardless of thread count. |
| `dipolar_geometry` | `true` | Generate random 2D dipolar couplings per realisation. `false` leaves the system uncoupled. |
| `dipolar_scale` | `1.0` | Multiplies the bare 1/r³ coupling. |
| `field_disorder` | `0.0` | Width of the Gaussian offset added to each Hz. |

Couplings cannot be listed explicitly yet — see
[architecture.md](architecture.md#extension-points).

## `[measurement]`

| Key | Default | Meaning |
|---|---|---|
| `reference_axis` | `"z"` | Axis b of the companion state A = Σ σ_b, prepared at t = 0 and propagated alongside. |
| `observables` | required | List of observable names. |

Observable names:

| Name | Columns | Meaning |
|---|---|---|
| `corr:<axis>` | `_re`, `_im` | Two-time correlator against the reference state — the echo signal |
| `mag:<spin>:<axis>` | one | Magnetisation of one spin, **1-based** |
| `norm` | two | Norms of both states; drift is numerical error |
| `energy` | one | ⟨H⟩, conserved, an independent check |

Axes are `x`, `y`, `z`. Include `norm` in production runs: it is free and it is
the cheapest evidence a result is numerically sound.

## `[numerics]`

| Key | Default | Meaning |
|---|---|---|
| `epsilon` | `1e-7` | Chebyshev coefficient cutoff. Smaller is more accurate and slower; cost grows roughly logarithmically. |

## `[[step]]` — the sequence

An ordered list. Four kinds.

### `measure`

Records every observable at the current time.

```toml
[[step]]
kind = "measure"
```

### `evolve`

Advances time under the Hamiltonian.

| Key | Default | Meaning |
|---|---|---|
| `tau` | required | Step length. Negative runs time backwards. |
| `steps` | `1` | Repetitions of the step. Cheaper than repeating the block: the expansion is built once. |
| `mode` | `"real"` | `real` or `imaginary`. |
| `autonormalize` | `false` | Rescale to unit norm after each step. Leave off for real time — it hides error rather than fixing it. |
| `hx`, `hy`, `hz` | baseline | Override fields for this step only. Components not named keep their `[field]` value. |

This is also how pulses are done realistically: a strong transverse field for a
short time. The flip angle is `H·τ`, so `hx = 10700` with
`tau = 0.000314159265` gives 3.3615 rad = **1.07π** — a deliberate 7%
over-rotation, which is the pulse-error physics the archived runs studied.

### `pulse`

An ideal instantaneous rotation, `exp(-i(θ/2) n·σ)`.

```toml
[[step]]
kind = "pulse"
rotation = "pi:x"        # every spin
# rotation = "pi/2:y:3"  # spin 3 only, 1-based
```

Angles: `pi`, `pi/2`, `pi/4`. Axes: `x`, `y`, `z`.

Use `pulse` for the idealised comparison and `evolve` for the realistic one;
subtracting the two separates flip-angle error from the error caused by couplings
acting during a finite pulse.

### `repeat`

Runs a body of steps `count` times. **Nests arbitrarily** — this is what makes
composite sequences expressible, and what the Fortran original refused to do.

```toml
[[step]]
kind = "repeat"
count = 20
body = [
  { kind = "evolve", tau = 0.05 },
  { kind = "evolve", tau = 0.000314159265, hx = 10700.0 },
  { kind = "evolve", tau = 0.05 },
  { kind = "measure" },
]
```

Bodies are inline tables in an array, and a body entry may itself be a `repeat`.

## Worked example: XY-8

The X and Y pulses differ only in which field component is driven, and nesting
gives the block structure.

```toml
[system]
spins = 12
initial_state = "random"

[field]
hx = 0.0
hy = 0.0
hz = 0.0

[ensemble]
realizations = 180
seed = 20240726
dipolar_geometry = true
field_disorder = 100.0

[measurement]
reference_axis = "x"
observables = ["corr:x", "norm", "energy"]

[[step]]
kind = "measure"

[[step]]
kind = "repeat"
count = 8                       # 8 blocks of XY-8 = 64 pulses
body = [
  { kind = "repeat", count = 2, body = [
      { kind = "evolve", tau = 0.05 },
      { kind = "evolve", tau = 0.000314159265, hx = 10700.0 },   # X
      { kind = "evolve", tau = 0.05 },
      { kind = "evolve", tau = 0.000314159265, hy = 10700.0 },   # Y
  ]},
  { kind = "measure" },
]
```

## Sweeps

There is no sweep driver; use the shell.

```bash
for eps in 1.00 1.03 1.07 1.10; do
  hx=$(python3 -c "print(10000 * $eps)")
  sed "s/hx = 10700.0/hx = $hx/" configs/cpmg8z.toml > /tmp/e$eps.toml
  ./build/spinsim run /tmp/e$eps.toml --out results/flip_$eps.csv
done
```

Command-line flags that override the file: `--realizations N`, `--threads N`,
`--out FILE`.

## Output

CSV with a self-describing header: `time`, then three columns per observable
column — `_mean`, `_sd`, `_stderr`.

```
time,corr_x_re_mean,corr_x_re_sd,corr_x_re_stderr,...
0,9.953018468,0.546761408695,0.122259567732
```

For error bars use **`_stderr`**, the uncertainty of the mean. `_sd` is the
spread between realisations — a property of the disorder, not of your estimate.
The Fortran original printed neither correctly; see [legacy.md](legacy.md).

## Migrating a legacy `.dat` script

```bash
python3 tools/dat2toml.py ../CPMGZ/CPMG8Z/nmrCPMG8.dat > configs/cpmg8.toml
```

The converter handles all 130 archived scripts. It cannot recover three things,
because they lived in Fortran source rather than in any data file, and it marks
each with a `TODO`: the measured correlator, the realisation count, and the
disorder strength. Fill them in by hand.
