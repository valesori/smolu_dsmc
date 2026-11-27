import math

import pytest

from smolu_dsmc.kernels import (
    blocking_factor,
    cyclization_rate,
    de_gennes_kernel,
    entangled_kernel,
    rouse_kernel,
)


def test_rouse_kernel():
    assert rouse_kernel(1, 1) == 2.0
    assert rouse_kernel(2, 4) == pytest.approx(0.75)


@pytest.mark.parametrize("gamma", [7e-6, 1e-3, 0.5, math.pi / 2])
def test_blocking_factor_limits(gamma):
    assert blocking_factor(0.0, gamma) == pytest.approx(1.0, rel=1e-6)
    assert blocking_factor(1e12, gamma) < 1e-6


def test_blocking_factor_decreases():
    values = [blocking_factor(x, 7e-6) for x in (0.0, 0.5, 10.0, 1e3, 1e5)]
    assert all(a > b for a, b in zip(values, values[1:], strict=False))


def test_entangled_kernel_large_masses():
    # Masses whose squared minimum overflows int32 (the original used int32).
    m = 50_000
    value = entangled_kernel(m, m + 1, density=1.0, gamma=7e-6, c0_a3=1.0)
    expected = blocking_factor(0.5 * m * m, 7e-6) * rouse_kernel(m, m + 1)
    assert value == pytest.approx(expected)
    assert 0.0 < value < rouse_kernel(m, m + 1)


def test_ring_kernels():
    assert de_gennes_kernel(1, 1) == 4.0
    assert cyclization_rate(1, 0.5) == 0.5
    assert cyclization_rate(4, 0.5) == pytest.approx(0.5 / 8)
