"""Command-line interface: ``python -m smolu_dsmc {frag,rings} [options]``."""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from . import fragmentation, fragmentation_gillespie, rings
from .output import FragmentationWriter, RingsWriter


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
        help="write a mass histogram every this many MC steps (default: %(default)s)",
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
    with FragmentationWriter(params, args.outdir, suffix=suffix) as writer:
        result = simulate(
            params,
            on_event=writer.on_event,
            on_snapshot=writer.on_snapshot,
            snapshot_every=args.hist_every,
        )
    writer.write_final(result.time, result.population)


def run_rings(args: argparse.Namespace) -> None:
    params = rings.RingsParams(
        ntot=args.ntot,
        density=args.density,
        prefactor=args.prefactor,
        alpha=args.alpha,
        time_max=args.time_max,
        seed=args.seed,
    )
    with RingsWriter(params, args.outdir, verbose=args.verbose) as writer:
        rings.simulate(params, on_event=writer.on_event)


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.model == "frag":
            run_fragmentation(args)
        else:
            run_rings(args)
    except ValueError as err:
        parser.error(str(err))
