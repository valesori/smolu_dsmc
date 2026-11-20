import math
from collections import Counter

import numpy as np
import pytest

from smolu_dsmc import FragmentationParams
from smolu_dsmc.cli import main
from smolu_dsmc.fragmentation_gillespie import (
    GroupedPopulation,
    harmonic,
    sample_inverse_index,
    simulate_gillespie,
)
from smolu_dsmc.population import UniformStream


@pytest.mark.parametrize("n", [0, 1, 2, 63, 4095, 4096, 10_000, 123_457])
def test_harmonic(n):
    exact = math.fsum(1.0 / i for i in range(1, n + 1))
    assert harmonic(n) == pytest.approx(exact, rel=1e-13, abs=0)


@pytest.mark.parametrize("k", [1, 2, 7, 40])
def test_sample_inverse_index(k):
    uniform = UniformStream(np.random.default_rng(k))
    draws = 100_000
    counts = Counter(sample_inverse_index(k, uniform) for _ in range(draws))
    assert set(counts) == set(range(1, k + 1))
    h = harmonic(k)
    for j, c in counts.items():
        expected = draws / (j * h)
        assert abs(c - expected) < 5 * math.sqrt(expected)


def test_grouped_population_bookkeeping():
    uniform = UniformStream(np.random.default_rng(2), block_size=101)
    pop = GroupedPopulation([1] * 300, max_mass=300)
    for _ in range(5000):
        polymers = [p for p in range(pop.n) if pop.mass(p) > 1]
        if pop.n > 1 and (not polymers or uniform() < 0.5):
            a = uniform.index(pop.n)
            b = (a + 1 + uniform.index(pop.n - 1)) % pop.n
            pop.set_mass(b, pop.mass(a) + pop.mass(b))
            pop.remove(a)
        else:
            pos = polymers[uniform.index(len(polymers))]
            m1 = 1 + uniform.index(pop.mass(pos) - 1)
            m2 = pop.mass(pos) - m1
            pop.set_mass(pos, m1)
            pop.add(m2)
    assert pop.total_mass == 300
    for g, members in enumerate(pop.groups):
        for slot, pos in enumerate(members):
            assert pop.mass(pos).bit_length() - 1 == g
            assert pop._slot[pos] == slot
    assert sum(len(g) for g in pop.groups) == pop.n
    inv, harm = pop.inv_total, pop.harm_total
    pop.recompute_totals()
    assert inv == pytest.approx(pop.inv_total) and harm == pytest.approx(pop.harm_total)
    # Bounds hold for every molecule.
    for pos in range(pop.n):
        m = pop.mass(pos)
        g = m.bit_length() - 1
        assert 1 / m <= pop.inv_bound[g] and harmonic(m - 1) <= pop.harm_bound[g]


def test_gillespie_equilibrium_mean_length():
    params = FragmentationParams(ntot=4000, kd=0.1, time_max=60.0, seed=11)
    samples = []
    simulate_gillespie(
        params,
        on_snapshot=lambda step, t, pop: samples.append(params.ntot / pop.n) if t > 20 else None,
        snapshot_every=200,
    )
    assert len(samples) > 50
    assert np.mean(samples) == pytest.approx((1 + math.sqrt(1 + 80)) / 2, rel=0.03)


def test_gillespie_entangled_conserves_mass():
    result = simulate_gillespie(FragmentationParams(ntot=500, time_max=5.0, gamma=7e-6))
    assert result.population.total_mass == 500
    assert result.time >= 5.0


def test_cli_gillespie_method(tmp_path):
    main(
        [
            "frag",
            "--method",
            "gillespie",
            "--ntot",
            "200",
            "--time-max",
            "2",
            "--outdir",
            str(tmp_path),
        ]
    )
    tag = "ntot200_Kd1.00e-03_density1.00e+00_gillespie"
    assert (tmp_path / f"mav_vs_t_{tag}_tmax2.00e+00.dat").exists()
    assert (tmp_path / f"histograms_{tag}").is_dir()
