"""Command-line interface: ``python -m smolu_dsmc {frag,rings} [options]``."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from . import fragmentation, fragmentation_gillespie, rings
from .output import FragmentationWriter, RingsWriter
from .provenance import RunRecord
from .recording import log_time_grid


def _add_common(parser: argparse.ArgumentParser, ntot: int) -> None:
    parser.add_argument(
        "--ntot", type=int, default=ntot, help="number of monomers (default: %(default)s)"
    )
    parser.add_argument(
        "--density", type=float, default=1.0, help="monomer number density (default: %(default)s)"
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=1.0,
        help="share of the time increment given to annealing, in [0, 1] (default: %(default)s)",
    )
    parser.add_argument(
        "--time-max",
        type=float,
        default=1e5,
        help="stop at this waiting time (default: %(default)g)",
    )
    parser.add_argument("--seed", type=int, default=123, help="random seed (default: %(default)s)")
    parser.add_argument(
        "--outdir", default=".", help="directory for output files (default: current directory)"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="smolu-dsmc",
        description="Direct Simulation Monte Carlo of Smoluchowski equations.",
    )
    sub = parser.add_subparsers(dest="model", required=True)

    frag = sub.add_parser("frag", help="annealing and fragmentation (Ref. 1)")
    _add_common(frag, ntot=10_000)
    frag.add_argument(
        "--kd", type=float, default=1e-3, help="dissociation constant K_d (default: %(default)g)"
    )
    frag.add_argument(
        "--gamma",
        type=float,
        default=None,
        help="max. annealing angle in rad; enables the entanglement kernel (e.g. 7e-6)",
    )
    frag.add_argument(
        "--hist-every",
        type=int,
        default=10_000,
        help="write a mass histogram every this many MC steps; not used with --n-points "
        "(default: %(default)s)",
    )
    frag.add_argument(
        "--n-points",
        type=int,
        default=None,
        help="record the time series and the histograms at this many log-spaced times "
        "from --t-min to --time-max, instead of after every event and every "
        "--hist-every steps",
    )
    frag.add_argument(
        "--t-min",
        type=float,
        default=None,
        help="first time of the --n-points grid (default: 0.01 * kd / density, about 1%% of "
        "the time over which the mean length doubles at the start)",
    )
    frag.add_argument(
        "--method",
        choices=["paper", "gillespie"],
        default="paper",
        help="'paper': the DSMC algorithm of Ref. 1, Appendix B (used for the published "
        "results); 'gillespie': Gillespie's stochastic simulation algorithm with thinning, "
        "which samples the same dynamics faster but is NOT the algorithm of the paper "
        "(--alpha is ignored). Default: %(default)s",
    )

    ring = sub.add_parser("rings", help="irreversible condensation with ring formation (Ref. 2)")
    _add_common(ring, ntot=1_000)
    ring.add_argument(
        "--prefactor",
        type=float,
        default=0.5,
        help="cyclization rate prefactor (default: %(default)s)",
    )
    ring.add_argument(
        "-v", "--verbose", action="store_true", help="print chain and ring counts after each event"
    )
    return parser


def run_fragmentation(args: argparse.Namespace) -> None:
    params = fragmentation.FragmentationParams(
        ntot=args.ntot,
        kd=args.kd,
        density=args.density,
        alpha=args.alpha,
        time_max=args.time_max,
        gamma=args.gamma,
        seed=args.seed,
    )
    gillespie = args.method == "gillespie"
    simulate = fragmentation_gillespie.simulate_gillespie if gillespie else fragmentation.simulate
    suffix = "_gillespie" if gillespie else ""

    if args.n_points is not None:
        t_min = args.t_min if args.t_min is not None else 0.01 * params.kd / params.density
        grid = log_time_grid(t_min, params.time_max, args.n_points)
        recording = {"mode": "time_grid", "n_points": args.n_points, "t_min": t_min}
    else:
        grid = None
        recording = {"mode": "every_event", "hist_every_steps": args.hist_every}

    writer = FragmentationWriter(params, args.outdir, suffix=suffix)
    record = RunRecord(
        writer.run_record_path,
        model="fragmentation",
        method="gillespie" if gillespie else "dsmc (Ref. 1, Appendix B)",
        params=params,
        recording=recording,
        argv=args.command_line,
    )
    record.start()
    try:
        with writer:
            if grid is None:
                result = simulate(
                    params,
                    on_event=writer.on_event,
                    on_snapshot=writer.on_snapshot,
                    snapshot_every=args.hist_every,
                )
            else:
                result = simulate(params, record_times=grid, on_record=writer.on_record)
        final_path = writer.write_final(result.time, result.population)
    except BaseException as err:
        record.finish("failed", error=repr(err))
        raise
    record.finish(
        final_time=result.time,
        steps=result.steps,
        n_molecules=result.population.n,
        outputs=[p.name for p in (writer.timeseries_path, writer.hist_dir, final_path)],
    )


def run_rings(args: argparse.Namespace) -> None:
    params = rings.RingsParams(
        ntot=args.ntot,
        density=args.density,
        prefactor=args.prefactor,
        alpha=args.alpha,
        time_max=args.time_max,
        seed=args.seed,
    )
    writer = RingsWriter(params, args.outdir, verbose=args.verbose)
    record = RunRecord(
        writer.run_record_path,
        model="rings",
        method="dsmc",
        params=params,
        recording={"mode": "every_event"},
        argv=args.command_line,
    )
    record.start()
    try:
        with writer:
            result = rings.simulate(params, on_event=writer.on_event)
    except BaseException as err:
        record.finish("failed", error=repr(err))
        raise
    record.finish(
        final_time=result.time,
        steps=result.steps,
        n_chains=result.state.n_chains,
        n_rings=result.state.n_rings,
        outputs=[writer.mav_path.name, writer.counts_path.name],
    )


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    # Recorded in the run record: the actual command line, or the arguments
    # passed to main() when called from Python.
    args.command_line = list(sys.argv) if argv is None else ["smolu-dsmc", *argv]
    try:
        if args.model == "frag":
            run_fragmentation(args)
        else:
            run_rings(args)
    except ValueError as err:
        parser.error(str(err))
