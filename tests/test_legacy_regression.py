"""``--entanglement-x reduced`` reproduces the package at commit 706e7fe file for file.

706e7fe is the package before the entanglement variable was changed. This is
not a check against the 2022 code (``legacy_2022/``, see
``test_legacy_2022.py``). The fixture was recorded with
``scripts/record_baseline.py record`` at 706e7fe; see that script. Identical
output needs the same platform (the blocking factor uses the C math library),
so the check is skipped on another platform.
"""

import importlib.util
import json
import platform
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "data" / "baseline_706e7fe.json"

spec = importlib.util.spec_from_file_location(
    "record_baseline", ROOT / "scripts/record_baseline.py"
)
record_baseline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(record_baseline)

baseline = json.loads(FIXTURE.read_text())
env = baseline["environment"]
same_platform = (sys.platform, platform.machine()) == (env["platform"], env["machine"])


def test_fixture_is_from_baseline_commit():
    assert baseline["git"]["commit"].startswith("706e7fe")
    assert baseline["git"]["uncommitted_changes"] is False


@pytest.mark.skipif(not same_platform, reason=f"baseline recorded on {env}")
@pytest.mark.parametrize("name", sorted(baseline["cases"]))
def test_matches_baseline(name):
    case = baseline["cases"][name]
    assert record_baseline.run_case(case["args"]) == case["outputs"], (
        f"{name} differs from the baseline recorded at 706e7fe in {env}"
    )
