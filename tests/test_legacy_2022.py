"""The scripts in legacy_2022/ are the 2022 code of commit fee1e58, unchanged."""

import hashlib
from pathlib import Path

import pytest

LEGACY = Path(__file__).resolve().parents[1] / "legacy_2022"

# SHA-256 of `git cat-file blob fee1e58:<name>`.
FEE1E58 = {
    "smoluchowski_dsmc_frag.py": "2c667789d87218643520626025b38ec3ea93b8003c150e8d9d7ced8e22ffb3da",
    "smoluchowski_dsmc_frag_entanglement.py": (
        "4c77c0a95785e115e448cd90806fed08516d9fbb339cc2e55e2e9f250ba82564"
    ),
    "smoluchowski_dsmc_rings.py": (
        "3f0c54c0e86a1fbb21cf16b0f9b1576b5dd92eda508a02d06760dd9ce7b1216c"
    ),
}


@pytest.mark.parametrize("name", sorted(FEE1E58))
def test_legacy_script_unchanged(name):
    digest = hashlib.sha256((LEGACY / name).read_bytes()).hexdigest()
    assert digest == FEE1E58[name], f"legacy_2022/{name} differs from commit fee1e58"
