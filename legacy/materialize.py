#!/usr/bin/env python3
"""Verify the legacy source archive or reconstruct one historical experiment."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
LEGACY_ROOT = REPO_ROOT / "legacy"
EXPERIMENTS_ROOT = LEGACY_ROOT / "experiments"
MANIFEST_PATH = LEGACY_ROOT / "source-manifest.csv"
SOURCE_ROOTS = (LEGACY_ROOT / "reference", LEGACY_ROOT / "variants")
SOURCE_SUFFIXES = {".f", ".h", ".mak", ".py"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_manifest() -> list[dict[str, str]]:
    with MANIFEST_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    expected = {"original_path", "canonical_path", "sha256"}
    if not rows or set(rows[0]) != expected:
        raise ValueError(f"unexpected columns in {MANIFEST_PATH}")
    return rows


def verify_manifest(rows: list[dict[str, str]]) -> None:
    originals: set[str] = set()
    canonical_digests: dict[str, str] = {}

    for row in rows:
        original = row["original_path"]
        canonical = row["canonical_path"]
        expected_digest = row["sha256"]
        if original in originals:
            raise ValueError(f"duplicate historical path: {original}")
        originals.add(original)

        canonical_path = (REPO_ROOT / canonical).resolve()
        if not any(canonical_path.is_relative_to(root.resolve()) for root in SOURCE_ROOTS):
            raise ValueError(f"canonical path escapes the source archive: {canonical}")
        if not canonical_path.is_file():
            raise FileNotFoundError(f"missing canonical source: {canonical}")

        actual_digest = canonical_digests.setdefault(canonical, sha256(canonical_path))
        if actual_digest != expected_digest:
            raise ValueError(f"checksum mismatch for {canonical}")

    archived_sources = {
        path.relative_to(REPO_ROOT).as_posix()
        for root in SOURCE_ROOTS
        for path in root.rglob("*")
        if path.is_file() and path.suffix in SOURCE_SUFFIXES
    }
    unreferenced = archived_sources - set(canonical_digests)
    if unreferenced:
        paths = ", ".join(sorted(unreferenced))
        raise ValueError(f"unreferenced canonical sources: {paths}")

    digest_to_path: dict[str, str] = {}
    for canonical, digest in canonical_digests.items():
        previous = digest_to_path.setdefault(digest, canonical)
        if previous != canonical:
            raise ValueError(f"duplicate canonical content: {previous}, {canonical}")

    print(
        f"Manifest OK: {len(rows)} historical paths, "
        f"{len(canonical_digests)} unique source files"
    )


def tracked_files(experiment: Path) -> list[Path]:
    relative = experiment.relative_to(REPO_ROOT).as_posix()
    if (REPO_ROOT / ".git").exists():
        result = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "ls-files", "-z", "--", relative],
            check=True,
            stdout=subprocess.PIPE,
        )
        return [
            REPO_ROOT / item.decode("utf-8")
            for item in result.stdout.split(b"\0")
            if item
        ]
    return [path for path in experiment.rglob("*") if path.is_file()]


def materialize(experiment_arg: str, output_arg: str, rows: list[dict[str, str]]) -> None:
    experiment = (Path.cwd() / experiment_arg).resolve()
    try:
        relative = experiment.relative_to(EXPERIMENTS_ROOT.resolve())
    except ValueError as error:
        raise ValueError(f"experiment must be inside {EXPERIMENTS_ROOT}") from error
    if len(relative.parts) != 2 or not experiment.is_dir():
        raise ValueError("select an experiment directory, for example legacy/experiments/rabi/Rabi")

    output = (Path.cwd() / output_arg).resolve()
    if output.exists():
        raise FileExistsError(f"output path already exists: {output}")
    if output.is_relative_to(experiment):
        raise ValueError("output path must not be inside the experiment directory")

    output.mkdir(parents=True)
    for source in tracked_files(experiment):
        if not source.is_file():
            continue
        destination = output / source.relative_to(experiment)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)

    historical_prefix = relative.name
    source_rows = [
        row
        for row in rows
        if row["original_path"].startswith(f"{historical_prefix}/")
    ]
    if not source_rows:
        raise ValueError(f"no historical sources found for {historical_prefix}")

    for row in source_rows:
        historical_relative = Path(row["original_path"]).relative_to(historical_prefix)
        destination = output / historical_relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO_ROOT / row["canonical_path"], destination)

    print(
        f"Materialized {relative.as_posix()} at {output} "
        f"with {len(source_rows)} source files"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment", nargs="?", help="experiment directory under legacy/experiments")
    parser.add_argument("--out", help="new output directory for the reconstructed experiment")
    parser.add_argument("--verify", action="store_true", help="verify manifest paths and SHA-256 checksums")
    args = parser.parse_args()
    if not args.verify and not args.experiment:
        parser.error("provide an experiment or use --verify")
    if args.experiment and not args.out:
        parser.error("--out is required when materializing an experiment")
    if args.out and not args.experiment:
        parser.error("--out requires an experiment")
    return args


def main() -> int:
    args = parse_args()
    try:
        rows = load_manifest()
        if args.verify:
            verify_manifest(rows)
        if args.experiment:
            materialize(args.experiment, args.out, rows)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
