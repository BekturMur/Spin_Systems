#!/usr/bin/env python3
"""Convert a legacy .dat pulse-sequence script into a spinsim TOML config.

This is a migration tool, run once per experiment, so it lives in Python rather
than in the C++ binary: it has to be forgiving of input the real parser should
reject, and none of it belongs in the simulation's hot path.

Usage:
    tools/dat2toml.py legacy/experiments/cpmg/CPMGZ/CPMG8Z/nmrCPMG8.dat \
        > configs/cpmg8z.toml

The legacy format is two nested languages:

  * a script file listing `@included.dat` steps, with `@CYCLE n` ... `@ENDCYCLE`
    around a repeated block;
  * per-step data files with keyword blocks terminated by `END`.

Both are matched by substring in chparsfPDDG.f, so `index(curline,'END')` closes
a block on any line containing those three letters anywhere. This reader is
deliberately stricter and complains rather than guessing.

What it cannot recover, and leaves as TODO comments:

  * the correlator being measured, which lived in chebsdPDDGnmr.f rather than in
    any data file (line 110: `addHSx` in CPMGZ, `addHSz` in Rabi_10_SzSz);
  * the realisation count, a `parameter (Nit=180)` in chebNMR2D-v6.f;
  * the disorder strength, `HzScale` in the same file.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import sys
from dataclasses import dataclass, field


@dataclass
class StepData:
    """One legacy data file: a field table plus timing."""

    name: str
    nspins: int = 0
    hx: float = 0.0
    hy: float = 0.0
    hz: float = 0.0
    uniform_fields: bool = True
    tau: float = 0.0
    nsteps: int = 0
    jz_scale: float = 1.0
    initial_state: str | None = None
    warnings: list[str] = field(default_factory=list)


def parse_data_file(path: pathlib.Path) -> StepData:
    """Read one .dat step file."""
    step = StepData(name=path.name)
    lines = path.read_text().splitlines()
    i = 0
    spin_fields: list[tuple[float, float, float]] = []

    def value_after(index: int) -> str:
        if index + 1 >= len(lines):
            raise ValueError(f"{path.name}: block at line {index + 1} has no value")
        return lines[index + 1].strip()

    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("HEISENBERG INTERACTION"):
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("END"):
                parts = [p for p in re.split(r"[,\s]+", lines[i].strip()) if p]
                if len(parts) >= 5:
                    # i, j, Jx, Jy, Jz -- j sets the system size, Jz the scale
                    # the driver later multiplies the dipolar couplings by.
                    step.nspins = max(step.nspins, int(parts[0]), int(parts[1]))
                    step.jz_scale = float(parts[4])
                i += 1

        elif line.startswith("SPIN :"):
            values = [p for p in re.split(r"[,\s]+", value_after(i)) if p]
            if len(values) < 9:
                raise ValueError(f"{path.name}: short SPIN row at line {i + 2}")
            # Columns are Hx0, Hx1, Omega_x, Hy0, ... ; only the static parts
            # (Hx0, Hy0, Hz0) are ever read by chparsfPDDG.f:205-207.
            spin_fields.append((float(values[0]), float(values[3]), float(values[6])))
            index = int(line.split(":")[1].split()[0])
            step.nspins = max(step.nspins, index)
            i += 2

        elif line.startswith("NUMBER OF FIELD STEPS"):
            step.nsteps = int(value_after(i))
            i += 2

        elif line.startswith("TIME STEP"):
            step.tau = float(value_after(i))
            i += 2

        elif line.startswith("NUMBER OF INTERMEDIATE TIME STEPS"):
            # chparsfPDDG.f:371 multiplies tau by this rather than subdividing.
            step.tau *= int(value_after(i))
            i += 2

        elif line.startswith("INITIAL STATE"):
            for keyword, name in (
                ("SPIN UP", "up"),
                ("SPIN DOWN", "down"),
                ("RANDOM", "random"),
            ):
                if keyword in line:
                    step.initial_state = name
            if "GROUND STATE" in line:
                step.warnings.append(
                    "INITIAL STATE: GROUND STATE is not supported; using random"
                )
            i += 1

        elif "SPECIAL PROPAGATION" in line:
            step.warnings.append(
                f"{path.name}: SPECIAL PROPAGATION (ideal gate) needs a manual "
                "[[step]] of kind = \"pulse\""
            )
            i += 1

        elif "IMAGINARY TIME PROPAGATION" in line:
            step.warnings.append(f'{path.name}: add mode = "imaginary" by hand')
            i += 1

        else:
            i += 1

    if spin_fields:
        step.hx, step.hy, step.hz = spin_fields[0]
        step.uniform_fields = all(f == spin_fields[0] for f in spin_fields)
        if not step.uniform_fields:
            step.warnings.append(
                f"{path.name}: fields differ between spins; only spin 1 was "
                "carried over, write the arrays by hand"
            )
    return step


@dataclass
class Script:
    steps: list[StepData]
    cycle_start: int | None
    cycle_count: int
    cycle_end: int | None


def parse_script(path: pathlib.Path) -> Script:
    """Read the top-level @-script."""
    steps: list[StepData] = []
    cycle_start: int | None = None
    cycle_end: int | None = None
    cycle_count = 1

    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("@CYCLE"):
            if cycle_start is not None:
                raise ValueError("nested @CYCLE is not representable")
            cycle_start = len(steps)
            tail = line[len("@CYCLE"):].strip()
            cycle_count = int(tail) if tail else 1
        elif line.startswith("@ENDCYCLE"):
            cycle_end = len(steps)
        elif line.startswith("@"):
            included = (path.parent / line[1:].strip()).resolve()
            if not included.exists():
                raise FileNotFoundError(f"{line} referenced from {path.name}")
            steps.append(parse_data_file(included))
    return Script(steps, cycle_start, cycle_count, cycle_end)


def render_step(step: StepData, indent: str) -> str:
    """A step file becomes an evolve, unless it only sets up the system."""
    if step.nsteps == 0:
        return f"{indent}# {step.name}: no field steps, setup only"
    parts = [f'kind = "evolve"', f"tau = {step.tau!r}", f"steps = {step.nsteps}"]
    for name, value in (("hx", step.hx), ("hy", step.hy), ("hz", step.hz)):
        if value != 0.0:
            parts.append(f"{name} = {value!r}")
    return f"{indent}{{ {', '.join(parts)} }},  # {step.name}"


def render(script: Script, source: pathlib.Path) -> str:
    if not script.steps:
        raise ValueError("the script includes no data files")

    nspins = max(s.nspins for s in script.steps)
    jz_scale = next((s.jz_scale for s in script.steps if s.jz_scale), 1.0)
    initial = next((s.initial_state for s in script.steps if s.initial_state), "random")

    warnings: list[str] = []
    for step in script.steps:
        warnings.extend(step.warnings)

    out: list[str] = []
    out.append(f"# Converted from {source} by tools/dat2toml.py.")
    out.append("#")
    out.append("# Check the TODOs below: they mark settings that the legacy code kept")
    out.append("# in Fortran source rather than in any data file, so no converter can")
    out.append("# recover them.")
    if warnings:
        out.append("#")
        for w in warnings:
            out.append(f"# WARNING: {w}")
    out.append("")
    out.append("[system]")
    out.append(f"spins = {nspins}")
    out.append(f'initial_state = "{initial}"')
    out.append("")

    setup = next((s for s in script.steps if s.nsteps == 0), None)
    out.append("[field]")
    out.append(f"hx = {setup.hx if setup else 0.0!r}")
    out.append(f"hy = {setup.hy if setup else 0.0!r}")
    out.append(f"hz = {setup.hz if setup else 0.0!r}")
    out.append("")

    out.append("[ensemble]")
    out.append("# TODO: legacy Nit, a `parameter` in chebNMR2D-v6.f:17")
    out.append("realizations = 180")
    out.append("seed = 1")
    out.append("dipolar_geometry = true")
    out.append(f"dipolar_scale = {jz_scale!r}")
    out.append("# TODO: legacy HzScale, chebNMR2D-v6.f:56")
    out.append("field_disorder = 100.0")
    out.append("")

    out.append("[measurement]")
    out.append("# TODO: the measured component lived in chebsdPDDGnmr.f:110, not here.")
    out.append('# CPMGZ used addHSx -> "x"; Rabi_10_SzSz used addHSz -> "z".')
    out.append('reference_axis = "z"')
    out.append('observables = ["corr:z", "norm"]')
    out.append("")

    out.append("[numerics]")
    out.append("epsilon = 1e-7")
    out.append("")

    out.append("[[step]]")
    out.append('kind = "measure"')
    out.append("")

    body = script.steps[script.cycle_start:script.cycle_end] if script.cycle_start is not None else []
    prologue = script.steps[:script.cycle_start] if script.cycle_start is not None else script.steps

    for step in prologue:
        rendered = render_step(step, "")
        if rendered.lstrip().startswith("#"):
            out.append(rendered)
            continue
        out.append("[[step]]")
        out.append(rendered.strip().rstrip(",").strip("{} ").replace(", ", "\n"))
        out.append("")

    if body:
        out.append("[[step]]")
        out.append('kind = "repeat"')
        out.append(f"count = {script.cycle_count}")
        out.append("body = [")
        for step in body:
            rendered = render_step(step, "  ")
            out.append(rendered)
        out.append('  { kind = "measure" },')
        out.append("]")
    out.append("")
    return "\n".join(out)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("script", type=pathlib.Path, help="legacy @-script .dat file")
    args = parser.parse_args()

    try:
        print(render(parse_script(args.script), args.script))
    except (ValueError, FileNotFoundError) as exc:
        print(f"dat2toml: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
