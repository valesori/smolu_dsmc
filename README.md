# smolu_dsmc

Direct Simulation Monte Carlo (DSMC) code in Python for Smoluchowski equations
with fragmentation, entanglement and ring formation. The models and the
algorithm are described in Refs. 1 and 2 below; the DSMC scheme follows
Appendix B of Ref. 1 (after Garcia et al., Physica A 143, 535 (1987)).

Three models are available:

- **Fragmentation** (Ref. 1): Smoluchowski coagulation equation with
  fragmentation, with annealing kernel `K_ij = (1/i + 1/j) / K_d` and
  fragmentation rate `F_ij = K_d K_ij / 2` (detailed balance).
- **Fragmentation with entanglement** (Ref. 1): the same, with the annealing
  and fragmentation rates multiplied by a blocking factor `g(x, gamma)` that
  accounts for the steric effect of entanglements.
- **Rings** (Ref. 2): irreversible coagulation (no fragmentation) with loop
  formation.

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

It was checked against the paper's algorithm (12–48 seeds, with and without
entanglement): the mean length agrees at all times within statistical error,
and the time-averaged equilibrium length is 14.671 ± 0.009 (Gillespie) vs
14.664 ± 0.010 (DSMC) for `ntot = 3000`, `K_d = 1e-2`, against 14.651 from
Eq. 14 of Ref. 1. The default remains `--method paper`. The algorithm is
described in the docstring of `smolu_dsmc/fragmentation_gillespie.py`.

### Output files

Fragmentation (`<tag>` = `ntot…_Kd…_density…[_gamma…][_gillespie]`):

- `mav_vs_t_<tag>_tmax….dat`: after each event, the waiting time, the mean
  filament length and the extent of reaction `1 - n/ntot`.
- `histograms_<tag>/n<step>_t<time>`: length histograms (length, count), every
  `--hist-every` Monte Carlo steps.
- `masses_final_<tag>_t….dat`: final lengths, padded with zeros to `ntot` entries.

Rings (`<tag>` = `ntot…_density…_pref…`):

- `mav_vs_t_<tag>.dat`: waiting time, mean mass of all molecules, of chains and of rings.
- `nchain_rings_<tag>.dat`: waiting time, number of chains and number of rings.

## Code layout

| Module | Contents |
| --- | --- |
| `smolu_dsmc/kernels.py` | Annealing, fragmentation, entanglement and cyclization rates |
| `smolu_dsmc/population.py` | O(1) bookkeeping of the molecules and random selection |
| `smolu_dsmc/fragmentation.py` | DSMC loop for annealing and fragmentation (Ref. 1, Appendix B) |
| `smolu_dsmc/fragmentation_gillespie.py` | Gillespie alternative, faster, not the algorithm of the paper |
| `smolu_dsmc/rings.py` | DSMC loop for annealing and ring formation |
| `smolu_dsmc/output.py` | Output files |
| `smolu_dsmc/cli.py` | Command-line interface |

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
