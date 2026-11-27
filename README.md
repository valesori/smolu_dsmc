# smolu_dsmc

Direct Simulation Monte Carlo (DSMC) code in Python for Smoluchowski equations
with fragmentation, entanglement and ring formation. The models and the
algorithm are described in Refs. 1 and 2 below; the DSMC scheme follows
Appendix B of Ref. 1, after Garcia et al. (Ref. 3). For
the fragmentation models, a faster Gillespie algorithm, which is not the one
used in Ref. 1, is also available (see below).

Three models are available:

- **Fragmentation** (Ref. 1): Smoluchowski coagulation equation with
  fragmentation, with annealing kernel `K_ij = (1/i + 1/j) / K_d` and
  fragmentation rate `F_ij = K_d K_ij / 2` (detailed balance).
- **Fragmentation with entanglement** (Ref. 1): the same, with the annealing
  and fragmentation rates multiplied by a blocking factor `g(x, gamma)` that
  accounts for the steric effect of entanglements.
- **Rings** (Ref. 2): irreversible coagulation (no fragmentation) with loop
  formation.

## Which code to use

| | Where | Use it for |
| --- | --- | --- |
| **2022 code** (legacy) | [`legacy_2022/`](legacy_2022/README.md) | Reproducing the published simulations of Ref. 1 exactly. The scripts of 2022, unchanged; they need NumPy < 2 (pinned environment in that folder). |
| **DSMC** | `smolu-dsmc frag` (`--method dsmc`, the default) | New runs. The same algorithm as Ref. 1, Appendix B, reimplemented: about 10 times faster, tested, with options for parameters and output. With entanglement, it uses the physical x by default (see "Entanglement"). |
| **Gillespie** | `smolu-dsmc frag --method gillespie` | New runs, faster still. A different algorithm that samples the same dynamics (see below); not used in Ref. 1. |

The package (DSMC and Gillespie) agrees with the 2022 code statistically, not
trajectory by trajectory: it draws random numbers differently.

## Units

The code is dimensionless and not tied to a particular system. It uses a unit
of time τ and a reference number density c₀. Masses (lengths) are numbers of
monomers (ULFs for Ref. 1).

- `--density` is the reduced monomer number density c̃ = c/c₀, where c counts
  all monomers per unit volume, including those in filaments. For Ref. 1,
  c̃ = 1 is the lowest experimental density, not one monomer per m³;
- the system volume is V = `ntot / density`, in units of 1/c₀, and two
  molecules anneal at rate K_ij/V;
- `--kd` (K_d) is in units of c₀;
- times are in units of τ, and the rates K_ij c₀, F_ij and the cyclization
  rate in units of 1/τ; `--gamma` is in radians;
- the entanglement variable x is the only place where the physical density
  enters (next section).

The physical values of τ and c₀ are fixed by the model being studied:

- **Ref. 1** sets τ = 1/F₁,₁, so that K_ij = (1/i + 1/j)/K_d and
  F_ij = (1/i + 1/j)/2 (Ref. 1, Eqs. 15–16). c₀ is the lowest concentration
  used (vimentin at 0.01 mg/mL), and K_d = 10⁻³ corresponds to the fitted
  1.0 × 10⁻⁵ mg/mL, or 6.1 pM (Ref. 1, Sec. II.B); hence c₀ = 6.1 nM of ULFs,
  or 3.67 × 10¹⁸ m⁻³. The early growth time scale is K_d/c: the mean length
  reaches about 2 at t ≈ K_d/c.
- **Ref. 2** (rings): the kernels are those of `smolu_dsmc/kernels.py`;
  density enters only the annealing rate, since cyclization involves a
  single chain. The units are those defined in Ref. 2 (not checked against
  that paper here).

## Entanglement

The blocking factor g(x, γ) takes x = c a³ l_min²/2 (Ref. 1, Eq. 4 and
Appendix D), where c is the physical density, l_min the number of ULFs of the
shorter filament and a the **axial length per ULF** (not the ULF diameter).
The 2022 code used for Ref. 1 computed x from the reduced density c̃ = c/c₀
instead, omitting the factor c₀a³; its fitted angle, γ = 7 × 10⁻⁶ rad, is
therefore not the physical angle. `--entanglement-x` selects how x is
computed:

| | `physical` (default) | `reduced` (legacy convention) |
| --- | --- | --- |
| x | c₀a³ · c̃ l_min²/2 = c a³ l_min²/2 | c̃ l_min²/2, as in the 2022 code |
| `--gamma` | physical maximum misalignment angle (rad) | angle fitted with the legacy convention (rad) |
| Ref. 1 fit | γ ≈ 8.82 × 10⁻³ rad (0.51°) | γ = 7 × 10⁻⁶ rad |
| c₀, a | `--reference-density` (m⁻³), `--axial-length` (m) | not used |
| Output files | `…_gamma…` | `…_gamma…_legacy` |

c₀a³ is computed from c₀ and a, which default to the values of Ref. 1:
c₀ = 6.1 nM = 3.67 × 10¹⁸ m⁻³ and a = 60 nm (nominal ULF length), giving
c₀a³ ≈ 7.93 × 10⁻⁴. The code does not change them; for another axial
increment (e.g. 49 nm) pass `--axial-length` explicitly.

Only x changes. The density used by the rest of the algorithm (volume,
annealing rates, `--t-min`), the fitted K_d and the blocking function g are
the same in both, and runs without `--gamma` are unaffected. `--gamma` is
never rescaled: the same value means a different angle in each. The run
record (`run_<tag>….json`) has an `entanglement` entry with the choice of x
and the angle and, for the physical x, c₀, a and c₀a³.

`--entanglement-x reduced` is the 2022 model in the new code: it gives the
same statistics as `legacy_2022/` but not the same trajectories (for those,
run `legacy_2022/`). It is useful to reuse the fitted angle with the faster
code or the Gillespie algorithm, and to compare the two conventions.

### Converting the fitted angle

The conversion of x is exact. That of the angle,

γ_phys ≃ γ_fit / (c₀a³) = 7 × 10⁻⁶ / 7.93 × 10⁻⁴ ≈ 8.82 × 10⁻³ rad,

is a small-angle approximation: for small γ, g depends on x and γ almost only
through γx. With these values g differs by at most 6 × 10⁻⁷ between the two
(`tests/test_entanglement_x.py`). Without `--entanglement-x reduced`,
`--gamma 7e-6` means a physical angle of 7 × 10⁻⁶ rad, for which
entanglements block almost nothing.

Examples (Ref. 1 parameters):

```bash
# Legacy convention, published fitted angle (output files: ..._legacy)
python -m smolu_dsmc frag --kd 1e-3 --gamma 7e-6 --entanglement-x reduced
# Physical x (default), corresponding physical angle (c0 and a: Ref. 1 defaults)
python -m smolu_dsmc frag --kd 1e-3 --gamma 8.82e-3
```

### Checks

- `legacy_2022/` holds the 2022 scripts byte for byte
  (`tests/test_legacy_2022.py`).
- `scripts/compare_legacy_2022.py` runs the 2022 scripts and the package
  side by side (in the environment of `legacy_2022/`). With ntot = 5000,
  K_d = 10⁻³, density 30 and 16 seeds, the mean lengths of the 2022 code and
  of the package (free; reduced x with γ = 7 × 10⁻⁶; physical x with
  γ = 8.82 × 10⁻³) agree within statistical error at 16 times from 10⁻³ to
  100 (largest |z| = 2.0 with entanglement, 2.5 without, over 48
  comparisons). The test is sensitive to entanglement: in the 2022 code it
  slows growth by 3–5.5 standard errors between t ≈ 0.05 and 5.
- `--entanglement-x reduced` gives, file for file, the output of the package
  before the physical x was introduced (commit 706e7fe;
  `tests/test_legacy_regression.py`).

## Installation

Requires Python >= 3.10 and NumPy >= 1.24.

```bash
pip install -e .            # or: pip install -e ".[test]" to also install pytest
```

## Usage

```bash
python -m smolu_dsmc frag  [--ntot 10000] [--kd 1e-3] [--density 1] [--alpha 1] \
                           [--time-max 1e5] [--hist-every 10000] [--seed 123] [--outdir .]
python -m smolu_dsmc frag --gamma 8.82e-3 [...]  # with entanglement (see above)
python -m smolu_dsmc frag --method gillespie [...]   # faster, not the paper's algorithm
python -m smolu_dsmc frag --n-points 200 [--t-min 1e-5] [...]   # record on a time grid
python -m smolu_dsmc rings [--ntot 1000] [--prefactor 0.5] [--density 1] [--alpha 1] \
                           [--time-max 1e5] [--seed 123] [--verbose] [--outdir .]
```

After `pip install`, the same commands are also available as `smolu-dsmc ...`.
The defaults are the parameter values of the original scripts, which are in
`legacy_2022/`.

The simulation can also be driven from Python:

```python
from smolu_dsmc import FragmentationParams
from smolu_dsmc.fragmentation import simulate

result = simulate(
    FragmentationParams(ntot=2000, kd=1e-2, time_max=50), on_event=lambda t, pop: print(t, pop.n)
)
```

### Gillespie algorithm for fragmentation (not the one used in the paper)

`--method gillespie` (Python:
`smolu_dsmc.fragmentation_gillespie.simulate_gillespie`) replaces the DSMC of
Ref. 1, Appendix B, with Gillespie's stochastic simulation algorithm (Refs. 4
and 5) combined with thinning (Ref. 6) and composition-rejection sampling
(Ref. 7). It samples the same dynamics, but it is **not the algorithm used for
the published results**. The DSMC rejects most of its attempts once filaments
are long (over 95% at equilibrium); the Gillespie version proposes reactions
with probabilities close to their rates, using the form of the kernels, and
draws exponential waiting times (`--alpha` is ignored). With `ntot = 1e5` it
is about 6 times faster up to `t = 20` and 14 times faster up to `t = 200`,
and the gain grows with filament length. Its output files carry a
`_gillespie` suffix.

It was checked against the DSMC with
`scripts/validate_gillespie.py` (about 20 s on 8 cores), which compares the
mean length at times from 0.01 to 200 over 24 seeds, with and without
entanglement (reduced x), and the time-averaged equilibrium length over 48 seeds. The two
agree within statistical error; the equilibrium length is 14.671 ± 0.009
(Gillespie) vs 14.664 ± 0.010 (DSMC) for `ntot = 3000`, `K_d = 1e-2`, against
14.651 from Eq. 14 of Ref. 1. The default remains `--method dsmc`. The
algorithm is described in the docstring of
`smolu_dsmc/fragmentation_gillespie.py`.

### Recording on a time grid

By default, the fragmentation model writes a line of the time series after
every accepted event, and a histogram every `--hist-every` Monte Carlo steps.
Long runs then produce very large files (the number of events grows with
`ntot` and with time: about 2 × 10⁶ lines up to t = 200 for `ntot = 1e5`), and
histograms taken every so many steps fall at irregular times, which differ
between the DSMC and the Gillespie algorithm.

`--n-points N` instead records the time series and the histograms at `N`
log-spaced times from `--t-min` to `--time-max`, with the state the system had
at each of these times. `--t-min` defaults to `0.01 * kd / density`, about 1%
of the time over which the mean length first doubles.

### Run records

Each run writes `run_<tag>….json` next to its output files, with the model,
the algorithm, all parameters (including the seed), the recording settings,
how the entanglement variable is computed, the package version, the git commit of the code and whether it had
uncommitted changes, the Python and NumPy versions, the command line, the
start and end times and the list of output files. Its `status` is `running`
while the run is in progress, then `finished` or `failed`.

### Output files

Fragmentation (`<tag>` = `ntot…_Kd…_density…[_gamma…[_legacy]][_gillespie]`):

- `mav_vs_t_<tag>_tmax….dat`: the waiting time, the mean filament length and
  the extent of reaction `1 - n/ntot`, after each event or at each grid time
  (`--n-points`).
- `histograms_<tag>/n<index>_t<time>`: length histograms (length, count),
  every `--hist-every` Monte Carlo steps (`<index>` is the step) or at each
  grid time (`<index>` is its position in the grid).
- `masses_final_<tag>_t….dat`: final lengths, padded with zeros to `ntot` entries.
- `run_<tag>_tmax….json`: the run record.

Rings (`<tag>` = `ntot…_density…_pref…`):

- `mav_vs_t_<tag>.dat`: waiting time, mean mass of all molecules, of chains and of rings.
- `nchain_rings_<tag>.dat`: waiting time, number of chains and number of rings.
- `run_<tag>.json`: the run record.

## Code layout

| Module | Contents |
| --- | --- |
| `smolu_dsmc/kernels.py` | Annealing, fragmentation, entanglement and cyclization rates |
| `smolu_dsmc/population.py` | O(1) bookkeeping of the molecules and random selection |
| `smolu_dsmc/fragmentation.py` | DSMC loop for annealing and fragmentation (Ref. 1, Appendix B) |
| `smolu_dsmc/fragmentation_gillespie.py` | Gillespie alternative, faster, not the algorithm of the paper |
| `smolu_dsmc/rings.py` | DSMC loop for annealing and ring formation |
| `smolu_dsmc/recording.py` | Recording the state at prescribed times |
| `smolu_dsmc/output.py` | Output files |
| `smolu_dsmc/provenance.py` | Run records (JSON) |
| `smolu_dsmc/cli.py` | Command-line interface |
| `scripts/validate_gillespie.py` | Comparison of the Gillespie and DSMC algorithms |
| `scripts/record_baseline.py` | Fixed-seed reference runs of the package (706e7fe) for `tests/test_legacy_regression.py` |
| `scripts/compare_legacy_2022.py` | Statistical comparison of the package with `legacy_2022/` |
| `legacy_2022/` | The 2022 scripts, unchanged, with their environment |

The tests include a check of the equilibrium mean length against the
analytical result, Eq. 14 of Ref. 1.

## Development

```bash
pip install -e ".[test]" pre-commit
pre-commit install          # runs ruff (lint and format) on each commit
pytest
```

See [CHANGELOG.md](CHANGELOG.md) for the changes from the original 2022 scripts.

## References

(1) “Fragmentation and Entanglement Limit Vimentin Intermediate Filament Assembly”, Quang D. Tran, Valerio Sorichetti, Gerard Pehau-Arnaudet, Martin Lenz, Cécile Leduc, Phys. Rev. X 13, 011014 (2023). Published version: https://doi.org/10.1103/PhysRevX.13.011014. Preprint (bioRxiv): https://doi.org/10.1101/2022.03.19.484978

(2) “Runaway Transition in Irreversible Polymer Condensation with Cyclization”, Maria Panoukidou, Simon Weir, Valerio Sorichetti, Yair Gutierrez Fosado, Martin Lenz, Davide Michieletto, Phys. Rev. Research 6, 023189 (2024). Published version: https://doi.org/10.1103/PhysRevResearch.6.023189. Preprint (arXiv): https://doi.org/10.48550/arXiv.2210.14010

(3) A. L. Garcia, C. van den Broeck, M. Aertsens and R. Serneels, “A Monte Carlo simulation of coagulation”, Physica A 143, 535–546 (1987). https://doi.org/10.1016/0378-4371(87)90164-6

(4) D. T. Gillespie, “Exact stochastic simulation of coupled chemical reactions”, J. Phys. Chem. 81, 2340–2361 (1977). https://doi.org/10.1021/j100540a008

(5) D. T. Gillespie, “An exact method for numerically simulating the stochastic coalescence process in a cloud”, J. Atmos. Sci. 32, 1977–1989 (1975).

(6) P. A. W. Lewis and G. S. Shedler, “Simulation of nonhomogeneous Poisson processes by thinning”, Naval Res. Logist. Q. 26, 403–413 (1979). https://doi.org/10.1002/nav.3800260304

(7) A. Slepoy, A. P. Thompson and S. J. Plimpton, “A constant-time kinetic Monte Carlo algorithm for simulation of large biochemical reaction networks”, J. Chem. Phys. 128, 205101 (2008). https://doi.org/10.1063/1.2919546
