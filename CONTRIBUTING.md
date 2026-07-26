# Contributing

Thank you for improving spinsim. Bug reports, reproducibility checks,
documentation fixes, and focused code changes are welcome.

## Development setup

The project requires CMake 3.28 or newer, a C++23 compiler, and LAPACK. The
first configure also downloads pinned Catch2 and toml++ dependencies.

```bash
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DSPINSIM_WERROR=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
python3 legacy/materialize.py --verify
```

If Ninja is unavailable, omit `-G Ninja`. Add or update tests for behavioral
changes. Numerical changes should state the expected tolerance and, where
possible, compare against the independent dense reference implementation.

## Project boundaries

- Maintained C++ code belongs in `include/`, `src/`, `apps/`, and `tests/`.
- User-facing examples belong in `configs/`; document new configuration keys in
  `docs/configuration.md`.
- Treat `legacy/` as an archive. Do not edit a hash-named canonical source in
  place without updating `legacy/source-manifest.csv` and its SHA-256 checksum.
- Do not commit generated builds, simulation outputs, TeX products, or copies of
  third-party papers.
- Preserve the sign, basis-ordering, and spin-normalization conventions in
  `docs/physics.md`, or call out an intentional compatibility break explicitly.

## Pull requests

Keep a pull request focused, explain the scientific or software motivation,
and include the commands used to verify it. CI builds with warnings as errors,
runs the complete test suite, verifies the legacy manifest, and checks that all
archived top-level scripts still convert to valid TOML.
