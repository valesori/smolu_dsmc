"""Compare the Gillespie algorithm with the DSMC algorithm of Ref. 1.

Reproduces the validation reported in the README:

1. Kinetics: the mean filament length <L>(t), averaged over seeds, at a few
   times, for the free kernel and for the entangled kernel.
2. Equilibrium: the time average of <L> over [20, 200] for the free kernel,
   averaged over seeds, compared with Eq. 14 of Ref. 1.

Both algorithms use the same seeds; with the default settings the numbers in
the README are reproduced exactly. Runs are spread over ``--jobs`` processes.

Usage:
    python scripts/validate_gillespie.py [--jobs 8] [--seeds-kinetics 24] [--seeds-equilibrium 48]
"""

from __future__ import annotations

import argparse
import math
import time
from multiprocessing import Pool

import numpy as np

from smolu_dsmc import FragmentationParams
from smolu_dsmc.fragmentation import simulate
from smolu_dsmc.fragmentation_gillespie import simulate_gillespie

METHODS = {"dsmc": simulate, "gillespie": simulate_gillespie}

KINETICS_TIMES = [0.01, 0.1, 0.5, 2.0, 10.0, 50.0, 200.0]
KINETICS_CASES = {
    "free": dict(ntot=3000, kd=1e-2, density=1.0, time_max=200.0),
    "entangled": dict(ntot=3000, kd=1e-3, density=5.0, gamma=7e-6, time_max=200.0),
}
KINETICS_FIRST_SEED = 1

EQUILIBRIUM_CASE = dict(ntot=3000, kd=1e-2, density=1.0, time_max=200.0)
EQUILIBRIUM_WINDOW = (20.0, 200.0)
EQUILIBRIUM_FIRST_SEED = 201


def run_kinetics(job: tuple[str, str, int]) -> tuple[str, str, list[float]]:
    case, method, seed = job
    params = FragmentationParams(seed=seed, **KINETICS_CASES[case])
    lengths: list[float] = []
    METHODS[method](
        params,
        record_times=KINETICS_TIMES,
        on_record=lambda i, t, pop: lengths.append(params.ntot / pop.n),
    )
    return case, method, lengths


def run_equilibrium(job: tuple[str, int]) -> tuple[str, float]:
    """Time average of <L> over EQUILIBRIUM_WINDOW for one run."""
    method, seed = job
    params = FragmentationParams(seed=seed, **EQUILIBRIUM_CASE)
    t0, t1 = EQUILIBRIUM_WINDOW
    history = [(0.0, 1.0)]
    METHODS[method](params, on_event=lambda t, pop: history.append((t, params.ntot / pop.n)))
    times, lengths = np.array(history).T
    # lengths[i] holds on [times[i], times[i+1]).
    edges = np.clip(np.append(times, t1), t0, t1)
    return method, float(np.sum(lengths * np.diff(edges)) / (t1 - t0))


def mean_and_error(values) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    return values.mean(axis=0), values.std(axis=0, ddof=1) / math.sqrt(len(values))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--jobs", type=int, default=8, help="parallel processes (default: 8)")
    parser.add_argument("--seeds-kinetics", type=int, default=24, help="seeds per case and method")
    parser.add_argument("--seeds-equilibrium", type=int, default=48, help="seeds per method")
    args = parser.parse_args()
    start = time.perf_counter()

    kinetics_jobs = [
        (case, method, seed)
        for case in KINETICS_CASES
        for method in METHODS
        for seed in range(KINETICS_FIRST_SEED, KINETICS_FIRST_SEED + args.seeds_kinetics)
    ]
    equilibrium_jobs = [
        (method, seed)
        for method in METHODS
        for seed in range(EQUILIBRIUM_FIRST_SEED, EQUILIBRIUM_FIRST_SEED + args.seeds_equilibrium)
    ]
    with Pool(args.jobs) as pool:
        kinetics = pool.map(run_kinetics, kinetics_jobs)
        equilibrium = pool.map(run_equilibrium, equilibrium_jobs)

    print(f"1. <L>(t), mean ± standard error over {args.seeds_kinetics} seeds")
    print(f"{'':20s}" + "".join(f"{f't = {t:g}':>16s}" for t in KINETICS_TIMES))
    for case, settings in KINETICS_CASES.items():
        print(f"   {case}: {settings}")
        for method in METHODS:
            values = [r for c, m, r in kinetics if c == case and m == method]
            mean, err = mean_and_error(values)
            cells = "".join(f"{f'{a:.2f} ± {e:.2f}':>16s}" for a, e in zip(mean, err, strict=True))
            print(f"   {method:17s}{cells}")

    t0, t1 = EQUILIBRIUM_WINDOW
    c, kd = EQUILIBRIUM_CASE["density"], EQUILIBRIUM_CASE["kd"]
    print(f"\n2. Time-averaged <L> over [{t0:g}, {t1:g}], {args.seeds_equilibrium} seeds")
    print(f"   {EQUILIBRIUM_CASE}")
    for method in METHODS:
        mean, err = mean_and_error([v for m, v in equilibrium if m == method])
        print(f"   {method:10s} {mean:.3f} ± {err:.3f}")
    print(f"   Eq. 14     {(1 + math.sqrt(1 + 8 * c / kd)) / 2:.3f}")
    print(f"\n({time.perf_counter() - start:.0f} s)")


if __name__ == "__main__":
    main()
