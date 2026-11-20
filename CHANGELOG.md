# Changelog

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
