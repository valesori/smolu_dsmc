"""Record, or check, fixed-seed reference runs of the fragmentation model.

    python scripts/record_baseline.py record tests/data/baseline_706e7fe.json
    python scripts/record_baseline.py check  tests/data/baseline_706e7fe.json

``record`` runs each case of :data:`CASES` through the command-line interface
and stores the SHA-256 of every output file (the run record excepted), with
the commit and environment the runs were made in. ``check`` reruns the cases
and compares. ``tests/test_legacy_regression.py`` does the same check.

The fixture in ``tests/data/baseline_706e7fe.json`` was recorded at commit
706e7fe, the package as it was before the entanglement variable was changed.
It checks that ``--entanglement-x reduced`` (the legacy convention) still
gives the package's earlier output, file for file. It says nothing about the
2022 code of ``legacy_2022/``, which the package reproduces only statistically
(different random number generator and bookkeeping).

Entangled cases are run with ``--entanglement-x reduced`` when the
command-line interface has that option (at 706e7fe the reduced x was the only
one). Their output files then carry a ``_legacy`` suffix, which is removed
before comparing names.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from smolu_dsmc.cli import build_parser, main  # noqa: E402
from smolu_dsmc.provenance import git_revision  # noqa: E402

CASES: dict[str, list[str]] = {
    "dsmc_entangled": "--ntot 3000 --gamma 7e-6 --time-max 200 --n-points 40 --seed 123",
    "dsmc_entangled_density5": (
        "--ntot 3000 --density 5 --gamma 7e-6 --time-max 200 --n-points 40 --seed 7"
    ),
    "dsmc_entangled_every_event": (
        "--ntot 500 --gamma 7e-6 --time-max 20 --hist-every 2000 --seed 11"
    ),
    "gillespie_entangled_density5": (
        "--method gillespie --ntot 3000 --density 5 --gamma 7e-6 --time-max 200 "
        "--n-points 40 --seed 123"
    ),
    "dsmc_free": "--ntot 3000 --kd 1e-2 --time-max 200 --n-points 40 --seed 123",
    "gillespie_free": (
        "--method gillespie --ntot 3000 --kd 1e-2 --time-max 200 --n-points 40 --seed 123"
    ),
}
CASES = {name: args.split() for name, args in CASES.items()}


def has_entanglement_x_option() -> bool:
    frag = build_parser()._subparsers._group_actions[0].choices["frag"]
    return any("--entanglement-x" in a.option_strings for a in frag._actions)


def case_argv(args: list[str]) -> list[str]:
    argv = ["frag", *args]
    if "--gamma" in args and has_entanglement_x_option():
        argv += ["--entanglement-x", "reduced"]
    return argv


def output_hashes(outdir: Path) -> dict[str, str]:
    """SHA-256 of every output file under ``outdir``, by relative path without ``_legacy``."""
    hashes = {}
    for path in sorted(outdir.rglob("*")):
        if path.is_file() and not path.name.startswith("run_"):
            name = str(path.relative_to(outdir)).replace("_legacy", "")
            hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


def run_case(args: list[str]) -> dict[str, str]:
    with tempfile.TemporaryDirectory() as tmp:
        main([*case_argv(args), "--outdir", tmp])
        return output_hashes(Path(tmp))


def environment() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "platform": sys.platform,
        "machine": platform.machine(),
    }


def record(path: Path) -> None:
    data = {
        "git": git_revision(),
        "environment": environment(),
        "cases": {name: {"args": args, "outputs": run_case(args)} for name, args in CASES.items()},
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def check(path: Path) -> bool:
    data = json.loads(path.read_text())
    ok = True
    for name, case in data["cases"].items():
        same = run_case(case["args"]) == case["outputs"]
        ok &= same
        print(f"{name}: {'identical' if same else 'DIFFERENT'}")
    return ok


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("action", choices=["record", "check"])
    parser.add_argument("fixture", type=Path)
    args = parser.parse_args()
    if args.action == "record":
        record(args.fixture)
    else:
        sys.exit(0 if check(args.fixture) else 1)
