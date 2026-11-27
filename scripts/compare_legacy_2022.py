"""Compare the 2022 fragmentation scripts (legacy_2022/) with the package, statistically.

    conda env create -f legacy_2022/environment.yml      # NumPy 1.26, see legacy_2022/README.md
    conda activate smolu-dsmc-legacy-2022
    python scripts/compare_legacy_2022.py [--seeds 16] [--jobs 8]

The 2022 scripts need NumPy < 2; the package also runs with NumPy 1.26, so one
environment runs both. For each seed, the 2022 script is copied into a
temporary directory with only its parameter lines changed (``ntot``, ``K_d``,
``density``, ``time_max``, ``myseed``, and ``nevery_hist`` to skip histograms),
and run. The package runs the same parameters with ``--method dsmc``: free, with
the legacy convention (``entanglement_x="reduced"``, gamma = 7e-6) and with the
physical x and the converted angle (gamma = 8.82e-3).

The mean length at log-spaced times is compared over seeds: the table gives
mean ± standard error and z = (difference) / (combined standard error). The
two codes use different random numbers, so they agree only statistically.
The parameters (density 30, K_d = 1e-3) are chosen so that entanglement
clearly slows growth between t ~ 0.05 and 5; the last column shows that
effect in the 2022 code itself.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from smolu_dsmc.fragmentation import FragmentationParams, simulate  # noqa: E402
from smolu_dsmc.recording import log_time_grid  # noqa: E402

LEGACY = ROOT / "legacy_2022"
PARAMS = dict(ntot=5000, kd=1e-3, density=30.0, time_max=100.0)
GRID = log_time_grid(1e-3, PARAMS["time_max"], 16)
PACKAGE_RUNS = {
    "free": {},
    "reduced": dict(gamma=7e-6, entanglement_x="reduced"),
    "physical": dict(gamma=8.82e-3),
}


def legacy_copy(source: str, seed: int) -> str:
    values = {
        "ntot": PARAMS["ntot"],
        "K_d": PARAMS["kd"],
        "density": PARAMS["density"],
        "time_max": PARAMS["time_max"],
        "myseed": seed,
        "nevery_hist": 10**12,
    }
    for name, value in values.items():
        source, n = re.subn(rf"^{name}=\S+", f"{name}={value!r}", source, flags=re.M)
        if n != 1:
            raise RuntimeError(f"parameter line {name!r} not found exactly once")
    return source


def run_legacy(job: tuple[str, int]) -> np.ndarray:
    script, seed = job
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "run.py"
        path.write_text(legacy_copy((LEGACY / script).read_text(), seed))
        subprocess.run([sys.executable, str(path)], cwd=tmp, check=True, capture_output=True)
        (series,) = Path(tmp).glob("mav_vs_t_*.dat")
        data = np.loadtxt(series, usecols=(0, 1), ndmin=2)
    last = np.searchsorted(data[:, 0], GRID, side="right") - 1  # last event before t
    return np.where(last >= 0, data[np.maximum(last, 0), 1], 1.0)


def run_package(job: tuple[str, int]) -> np.ndarray:
    kind, seed = job
    params = FragmentationParams(seed=seed, **PARAMS, **PACKAGE_RUNS[kind])
    lengths: list[float] = []
    simulate(
        params,
        record_times=GRID,
        on_record=lambda i, t, pop: lengths.append(pop.total_mass / pop.n),
    )
    return np.array(lengths)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--seeds", type=int, default=16)
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()
    if int(np.__version__.split(".")[0]) >= 2:
        sys.exit("The 2022 scripts need NumPy < 2: use legacy_2022/environment.yml.")

    seeds = range(args.seeds)
    with ProcessPoolExecutor(args.jobs) as ex:
        legacy = {
            name: ex.map(run_legacy, [(script, s) for s in seeds])
            for name, script in [
                ("2022 free", "smoluchowski_dsmc_frag.py"),
                ("2022 entangled", "smoluchowski_dsmc_frag_entanglement.py"),
            ]
        }
        package = {
            f"pkg {kind}": ex.map(run_package, [(kind, s) for s in seeds]) for kind in PACKAGE_RUNS
        }
        rows = {name: np.array(list(runs)) for name, runs in {**legacy, **package}.items()}

    mean = {k: v.mean(0) for k, v in rows.items()}
    err = {k: v.std(0, ddof=1) / np.sqrt(len(v)) for k, v in rows.items()}

    def z(a: str, b: str) -> np.ndarray:
        return (mean[a] - mean[b]) / np.hypot(err[a], err[b])

    print(f"{PARAMS}, {args.seeds} seeds; mean length ± standard error")
    tests = {
        "z(red.)": ("2022 entangled", "pkg reduced"),
        "z(phys.)": ("2022 entangled", "pkg physical"),
        "z(free)": ("2022 free", "pkg free"),
        "z(2022 ent.-free)": ("2022 entangled", "2022 free"),
    }
    print(f"{'t':>8s}" + "".join(f"{k:>17s}" for k in rows) + "".join(f"{k:>18s}" for k in tests))
    for i, t in enumerate(GRID):
        cells = "".join(f"{mean[k][i]:10.1f} ±{err[k][i]:5.1f}" for k in rows)
        zs = "".join(f"{z(*pair)[i]:18.2f}" for pair in tests.values())
        print(f"{t:8.3g}{cells}{zs}")
    for name, pair in list(tests.items())[:3]:
        print(f"max |{name}| = {np.abs(z(*pair)).max():.2f}")


if __name__ == "__main__":
    main()
