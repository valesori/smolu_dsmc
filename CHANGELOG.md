# Changelog

## [Unreleased]

### Changed (breaking)

- Entanglement: the variable x of the blocking factor is now computed from the
  physical density, x = c₀a³ c̃ l_min²/2 with a the axial length per ULF, and
  `--gamma` is the physical maximum misalignment angle. The 2022 code used
  x = c̃ l_min²/2. Without further options, `--gamma 7e-6` therefore no longer
  reproduces the fit of Ref. 1: use `--entanglement-x reduced`, or the
  corresponding physical angle `--gamma 8.82e-3` (README, "Entanglement").
  The density used elsewhere, K_d and the blocking function are unchanged.
  Files of runs made with `--gamma` between the 2026 rewrite and this change
  have no suffix but use the reduced x; their run record gives the commit.
- `--method paper` is now `--method dsmc` (still the default).

### Removed

- The top-level wrappers `smoluchowski_dsmc_frag.py`,
  `smoluchowski_dsmc_frag_entanglement.py` and `smoluchowski_dsmc_rings.py`.
  They had the names of the 2022 scripts but ran the package, with other
  inputs and random numbers. Use `smolu-dsmc frag` / `smolu-dsmc rings`
  (`python -m smolu_dsmc ...`), or the 2022 scripts in `legacy_2022/`.

### Added

- `legacy_2022/`: the three 2022 scripts (commit `fee1e58`) unchanged, with a
  pinned environment (Python 3.11, NumPy 1.26.4; they fail with NumPy >= 2)
  and a README. They reproduce the simulations of Ref. 1 exactly.
  `tests/test_legacy_2022.py` checks that they stay identical to `fee1e58`;
  `scripts/compare_legacy_2022.py` compares them with the package
  statistically.
- `--entanglement-x {physical,reduced}`, `--reference-density` and
  `--axial-length` (Python: `FragmentationParams.entanglement_x`,
  `reference_density`, `axial_length`); c₀ and a default to Ref. 1 (6.1 nM,
  60 nm). `reduced` is the legacy convention of the 2022 code; its output
  files carry a `_legacy` suffix. It reproduces the package at commit
  706e7fe file for file (`tests/test_legacy_regression.py`,
  `scripts/record_baseline.py`). Run records have an `entanglement` entry.
- `smolu_dsmc.fragmentation_gillespie` and `smolu-dsmc frag --method gillespie`:
  Gillespie's stochastic simulation algorithm with thinning and
  composition-rejection sampling (references in the README), which samples
  the same dynamics as the DSMC of Ref. 1 but is **not the algorithm used for
  the published results**. It avoids most rejected attempts and is 6–14 times
  faster for `ntot = 1e5`. Output files carry a `_gillespie` suffix; the
  default method is still the DSMC.
- `--n-points` and `--t-min` (fragmentation): record the time series and the
  histograms at log-spaced times instead of after every event, which keeps
  output files small and makes histograms comparable between runs and
  algorithms. Library: `record_times` and `on_record` arguments of both
  `simulate` functions (`smolu_dsmc.recording`).
- Run records: every run writes `run_<tag>….json` with its parameters, seed,
  algorithm, code version and git commit, environment and outputs
  (`smolu_dsmc.provenance`).
- `scripts/validate_gillespie.py`, which reproduces the comparison of the
  Gillespie and DSMC algorithms reported in the README.

### Changed

- README: the dimensionless units used by the code, how they map onto Refs. 1
  and 2, and the two entanglement conventions.
- Package docstring and description cover both algorithms and cite the
  published versions of Refs. 1 and 2.

## [1.0.0] - unreleased

Refactoring of the original 2022 scripts into the `smolu_dsmc` package. The
DSMC algorithm is unchanged (Appendix B of Ref. 1); its results agree with the
original scripts within statistical error.

### Changed

- The three scripts are now one package with separate modules for kernels,
  molecule bookkeeping, simulation loops, output and command-line interface.
  `python -m smolu_dsmc {frag,rings}` (or `smolu-dsmc`) replaces them, with the
  original parameter values as defaults. The old script names remain as
  wrappers.
- Requires Python >= 3.10 and NumPy >= 1.24; star imports are gone.
- Molecules are kept in a compact list instead of a zero-padded array scanned
  at every step. A run with `ntot = 10000` is about 10 times faster.
- Random numbers come from `numpy.random.default_rng`: runs with the same seed
  do not reproduce the old trajectories, only their statistics.
- The waiting-time series of the entanglement model has the same three
  columns as the plain fragmentation model (the extent of reaction was
  missing).
- Rings: the initial estimate of the maximum cyclization rate is its exact
  value (it was half of it), and printing counts after each event is now
  optional (`--verbose`).

### Fixed

- The entanglement script failed with NumPy >= 2, where `from numpy import *`
  shadows the built-in `min`.
- Masses were stored as `int32`, so `min(i, j)**2` in the entanglement kernel
  overflowed when both masses exceeded 46341.
- The blocking factor `g(x, gamma)` is evaluated without numerical
  cancellation: relative error ~1e-12 instead of ~1e-6 at `gamma = 7e-6`.
- Input checks: `K_d` must be > 0 (0 was accepted), `gamma` in (0, pi/2]
  (0 was accepted), `alpha` in [0, 1].
- A system reduced to a single filament no longer divides by zero.

### Added

- Tests (`pytest`), including the equilibrium mean length against Eq. 14 of
  Ref. 1.
- `pyproject.toml`, and a pre-commit hook running ruff.
