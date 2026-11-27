# legacy_2022: the original scripts

The three scripts of this folder are the code of December 2022, unchanged
(commit `fee1e58`, the last version before the 2026 rewrite;
`tests/test_legacy_2022.py` checks that they are identical to it):

| Script | Model | Reference |
| --- | --- | --- |
| `smoluchowski_dsmc_frag.py` | Annealing and fragmentation | Ref. 1 |
| `smoluchowski_dsmc_frag_entanglement.py` | The same, with entanglement | Ref. 1 |
| `smoluchowski_dsmc_rings.py` | Irreversible condensation with rings | Ref. 2 (cited in its header) |

The two fragmentation scripts are the code used for the simulations of Ref. 1
(Tran, Sorichetti et al., Phys. Rev. X 13, 011014 (2023)). Use them to
reproduce those simulations exactly. Whether the rings script produced the
simulations of Ref. 2 is not documented here. For new simulations, use the
package (`smolu-dsmc`, top-level [README](../README.md)), which is faster and
tested.

## Environment

The scripts need **NumPy < 2**. With NumPy 2, `from numpy import *` replaces
Python's built-in `min`, so `min(mi, mj)` becomes `numpy.min(mi, axis=mj)` and
the entanglement script fails at its first annealing. They have been run with
Python 3.11.18 and NumPy 1.26.4 (macOS arm64, 2026-11-27), which this folder
pins:

```bash
conda env create -f legacy_2022/environment.yml     # or: pip install -r legacy_2022/requirements.txt
conda activate smolu-dsmc-legacy-2022
```

The environment used in 2022 is not documented. With the same seed, a run
gives identical files when repeated. Other NumPy 1.x versions draw the same
random numbers (`numpy.random.seed` and `rand` have not changed); the last
digits of `exp`, `sin` and `cos` may differ between platforms, which can
change a trajectory, not its statistics.

## Running

Each script takes no options: the parameters are set in the lines at its top
(`K_d`, `ntot`, `density`, `gamma`, `time_max`, `myseed`, `Alpha`,
`nevery_hist`, and `glob_pref` for the rings). To run other parameters, copy
the script, edit those lines and run the copy in an empty directory, where it
writes its output:

```bash
mkdir run && cp legacy_2022/smoluchowski_dsmc_frag_entanglement.py run/
cd run    # edit the parameter lines of the copy, then:
python smoluchowski_dsmc_frag_entanglement.py
```

Fragmentation output (`<tag>` = `ntot…_Kd…_density…[_gamma…]`):

- `mav_vs_t_<tag>_tmax….dat`: waiting time and mean filament length after
  each event (plus the extent of reaction in `smoluchowski_dsmc_frag.py`);
- `histograms_<tag>/n<step>_t<time>`: length histograms every `nevery_hist`
  Monte Carlo steps;
- `masses_final_<tag>_t….dat`: final lengths, padded with zeros.

## Conventions

- **Units**: the same as the package. Densities and K_d are in units of c₀,
  the lowest experimental density of Ref. 1, and times in units of
  τ = 1/F₁,₁: K_ij = (1/i + 1/j)/K_d, F_ij = (1/i + 1/j)/2.
- **Entanglement**: x = (c/c₀) l_min²/2, from the reduced density, and
  `gamma = 7e-6` is the angle fitted with that x. The package calls this the
  legacy convention (`--entanglement-x reduced`); by default it computes x
  from the physical density, where the same fit corresponds to
  γ ≈ 8.82 × 10⁻³ rad (top-level README, "Entanglement").

## Relation to the package

`smolu-dsmc frag --method dsmc --entanglement-x reduced` implements the same
model with the same algorithm (Ref. 1, Appendix B): the same rates and
acceptance rule, and the same probabilities of choosing pairs and filaments.
It does not reproduce the trajectories of these scripts, because it draws
random numbers differently. `scripts/compare_legacy_2022.py` (run in this
environment) compares the two statistically: for ntot = 5000, K_d = 10⁻³,
density 30 and 16 seeds, the mean lengths agree within statistical error
from t = 10⁻³ to 100, with and without entanglement.

The package also differs where these scripts are inaccurate:

- masses are stored as `int32`, so `min(mi, mj)**2` overflows when both
  filaments exceed 46 341 ULFs (possible only for ntot > 92 682);
- g(x, γ) is evaluated as written, which loses about 6 significant digits to
  rounding at γ = 7 × 10⁻⁶ (relative error ~10⁻⁶);
- a system reduced to a single filament divides by zero.

## History

The scripts were first committed on 2022-12-21 (`af81ccb`) and revised the
same day up to `fee1e58`. Before `fee1e58` ("Fixed units"), the fragmentation
scripts used K_ij = (1/i + 1/j) and F_ij = K_d (1/i + 1/j)/2, i.e. a time unit
1/K_d times longer than the one above, and other default parameters.
