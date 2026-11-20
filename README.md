# smolu_dsmc

Direct Simulation Monte Carlo (DSMC) code in Python for Smoluchowski equations
with fragmentation, entanglement and ring formation. The models and the
algorithm are described in Refs. 1 and 2 below; the DSMC scheme follows
Appendix B of Ref. 1 (after Garcia et al., Physica A 143, 535 (1987)). For
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

## Units

The code is dimensionless and not tied to a particular system. It uses a unit
of time τ, a unit of number density c₀ and a unit of length b, related by
c₀b³ = 1 (one monomer per volume b³). Masses (lengths) are numbers of
monomers. Every quantity that involves density uses these same units:

- `--density` is the monomer number density in units of c₀, i.e. c/c₀ = c b³;
- the system volume is V = `ntot / density` (in units of b³), and two
  molecules anneal at rate K_ij/V;
- `--kd` (K_d) is in units of c₀;
- the entanglement variable x = c b³ l_min²/2 is `density * l_min**2 / 2`;
- times are in units of τ, and the rates K_ij c₀, F_ij and the cyclization
  rate in units of 1/τ; `gamma` is in radians.

The physical values of τ and c₀ (and hence b) are fixed by the model being
studied:

- **Ref. 1** sets τ = 1/F₁,₁, so that K_ij = (1/i + 1/j)/K_d and
  F_ij = (1/i + 1/j)/2 (Ref. 1, Eqs. 15–16). It expresses the experimental
  densities relative to the lowest concentration used (vimentin at
  0.01 mg/mL), and K_d = 10⁻³ corresponds to the fitted 1.0 × 10⁻⁵ mg/mL.
  The early growth time scale is K_d/c: the mean length reaches about 2 at
  t ≈ K_d/c.
- **Ref. 2** (rings): the kernels are those of `smolu_dsmc/kernels.py`;
  density enters only the annealing rate, since cyclization involves a
  single chain. The units are those defined in Ref. 2 (not checked against
  that paper here).

**Entanglement and the value of `gamma`.** The code computes the entanglement
variable from the dimensionless density c̃ = c/c₀ as x = c̃ l_min²/2. In
physical units x = c b³ l_min²/2 = c̃ (c₀b³) l_min²/2 (Ref. 1, Appendix D,
Eqs. 31–36), so the code omits the factor c₀b³; Ref. 1's Eq. 4, written
with dimensionless densities, omits it as well. This only changes how
`gamma` is read: for the small angles of interest the blocking factor
g(x, γ) depends only on the product γx (rescaling x by any factor and γ by its
inverse changes g by less than 3 × 10⁻⁷ for γ = 7 × 10⁻⁶), so the value used
in the code is γ_code = γ × (c₀b³), where γ is the bare angle.

Example for Ref. 1 (to be checked): K_d = 1.0 × 10⁻⁵ mg/mL is equivalent to
6.1 pM (Ref. 1, Sec. II.B), so c₀ = 0.01 mg/mL is 6.1 nM of ULFs, or
3.7 × 10¹⁸ m⁻³. Taking b as the ULF length, 60 nm (in Appendix D, b l_min is
a filament length), c₀b³ ≈ 7.9 × 10⁻⁴, and the fitted γ_code = 7 × 10⁻⁶ rad
corresponds to a bare angle γ ≈ 8.8 × 10⁻³ rad ≈ 0.5° (about 0.3–0.8° for
b = 55–65 nm and γ_code = (7 ± 2) × 10⁻⁶).

## Installation

Requires Python >= 3.10 and NumPy >= 1.24.

```bash
pip install -e .            # or: pip install -e ".[test]" to also install pytest
```

## Usage

```bash
python -m smolu_dsmc frag  [--ntot 10000] [--kd 1e-3] [--density 1] [--alpha 1] \
                           [--time-max 1e5] [--hist-every 10000] [--seed 123] [--outdir .]
python -m smolu_dsmc frag --gamma 7e-6 [...]     # with entanglement
python -m smolu_dsmc frag --method gillespie [...]   # faster, not the paper's algorithm
python -m smolu_dsmc frag --n-points 200 [--t-min 1e-5] [...]   # record on a time grid
python -m smolu_dsmc rings [--ntot 1000] [--prefactor 0.5] [--density 1] [--alpha 1] \
                           [--time-max 1e5] [--seed 123] [--verbose] [--outdir .]
```

After `pip install`, the same commands are also available as `smolu-dsmc ...`.
The defaults are the parameter values of the original scripts. The scripts
`smoluchowski_dsmc_frag.py`, `smoluchowski_dsmc_frag_entanglement.py` (which
sets `--gamma 7e-6`) and `smoluchowski_dsmc_rings.py` are kept as thin
wrappers and accept the same options.

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
Ref. 1, Appendix B, with Gillespie's stochastic simulation algorithm (Refs. 3
and 4) combined with thinning (Ref. 5) and composition-rejection sampling
(Ref. 6). It samples the same dynamics, but it is **not the algorithm used for
the published results**. The DSMC rejects most of its attempts once filaments
are long (over 95% at equilibrium); the Gillespie version proposes reactions
with probabilities close to their rates, using the form of the kernels, and
draws exponential waiting times (`--alpha` is ignored). With `ntot = 1e5` it
is about 6 times faster up to `t = 20` and 14 times faster up to `t = 200`,
and the gain grows with filament length. Its output files carry a
`_gillespie` suffix.

It was checked against the paper's algorithm with
`scripts/validate_gillespie.py` (about 20 s on 8 cores), which compares the
mean length at times from 0.01 to 200 over 24 seeds, with and without
entanglement, and the time-averaged equilibrium length over 48 seeds. The two
agree within statistical error; the equilibrium length is 14.671 ± 0.009
(Gillespie) vs 14.664 ± 0.010 (DSMC) for `ntot = 3000`, `K_d = 1e-2`, against
14.651 from Eq. 14 of Ref. 1. The default remains `--method paper`. The
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
the package version, the git commit of the code and whether it had
uncommitted changes, the Python and NumPy versions, the command line, the
start and end times and the list of output files. Its `status` is `running`
while the run is in progress, then `finished` or `failed`.

### Output files

Fragmentation (`<tag>` = `ntot…_Kd…_density…[_gamma…][_gillespie]`):

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

(3) D. T. Gillespie, “Exact stochastic simulation of coupled chemical reactions”, J. Phys. Chem. 81, 2340–2361 (1977). https://doi.org/10.1021/j100540a008

(4) D. T. Gillespie, “An exact method for numerically simulating the stochastic coalescence process in a cloud”, J. Atmos. Sci. 32, 1977–1989 (1975).

(5) P. A. W. Lewis and G. S. Shedler, “Simulation of nonhomogeneous Poisson processes by thinning”, Naval Res. Logist. Q. 26, 403–413 (1979). https://doi.org/10.1002/nav.3800260304

(6) A. Slepoy, A. P. Thompson and S. J. Plimpton, “A constant-time kinetic Monte Carlo algorithm for simulation of large biochemical reaction networks”, J. Chem. Phys. 128, 205101 (2008). https://doi.org/10.1063/1.2919546
