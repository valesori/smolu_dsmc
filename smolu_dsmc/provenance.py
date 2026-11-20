"""Run records: a JSON file describing how a set of output files was produced.

Each run started from the command line writes ``run_<tag>.json`` next to its
output files, with the model, the algorithm, all parameters (including the
seed), the recording settings, the package version, the git commit of the code
(and whether it had uncommitted changes), the Python and NumPy versions, the
command line, start and end times, and the list of output files. The file is
written when the run starts (``"status": "running"``) and updated when it
ends (``"finished"`` or ``"failed"``).
"""

from __future__ import annotations

import dataclasses
import datetime
import json
import platform
import subprocess
import sys
from importlib import metadata
from pathlib import Path
from typing import Any

import numpy as np

_PACKAGE_DIR = Path(__file__).resolve().parent


def package_version() -> str:
    try:
        return metadata.version("smolu-dsmc")
    except metadata.PackageNotFoundError:
        from . import __version__

        return __version__


def git_revision(path: Path = _PACKAGE_DIR) -> dict[str, Any] | None:
    """Commit and uncommitted-changes flag of the git checkout holding ``path``.

    None if the code is not in a git checkout (e.g. a regular pip install) or
    git is not available.
    """

    def git(*args: str) -> str:
        return subprocess.run(
            ["git", "-C", str(path), *args],
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
        ).stdout.strip()

    try:
        commit = git("rev-parse", "HEAD")
        dirty = bool(git("status", "--porcelain", "--untracked-files=no"))
    except (OSError, subprocess.SubprocessError):
        return None
    return {"commit": commit, "uncommitted_changes": dirty}


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


class RunRecord:
    """The JSON run record of one simulation."""

    def __init__(
        self,
        path: Path | str,
        *,
        model: str,
        method: str,
        params: Any,
        recording: dict[str, Any],
        argv: list[str] | None = None,
    ):
        self.path = Path(path)
        self.data: dict[str, Any] = {
            "status": "running",
            "model": model,
            "method": method,
            "parameters": dataclasses.asdict(params),
            "recording": recording,
            "code": {
                "package": "smolu-dsmc",
                "version": package_version(),
                "git": git_revision(),
            },
            "environment": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "platform": platform.platform(),
            },
            "command": argv if argv is not None else sys.argv,
            "started": _now(),
        }

    def start(self) -> None:
        self._write()

    def finish(self, status: str = "finished", **results: Any) -> None:
        self.data["status"] = status
        self.data["ended"] = _now()
        self.data["results"] = results
        self._write()

    def _write(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.data, indent=2, default=str) + "\n")
