import math

import numpy as np
import pytest

from smolu_dsmc import fragmentation, rings
from smolu_dsmc.cli import main


def mean_length_at_equilibrium(c, kd):
    """Ref. 1, Eq. 14."""
    return (1 + math.sqrt(1 + 8 * c / kd)) / 2


@pytest.mark.parametrize("alpha", [1.0, 0.5])
def test_equilibrium_mean_length(alpha):
    params = fragmentation.FragmentationParams(
        ntot=4000, kd=0.1, density=1.0, alpha=alpha, time_max=60.0, seed=7
    )
    samples = []
    fragmentation.simulate(
        params,
        on_snapshot=lambda step, t, pop: samples.append(params.ntot / pop.n) if t > 20 else None,
        snapshot_every=500,
    )
    assert len(samples) > 50
    assert np.mean(samples) == pytest.approx(mean_length_at_equilibrium(1.0, 0.1), rel=0.03)


def test_entangled_run_conserves_mass():
    params = fragmentation.FragmentationParams(ntot=500, time_max=5.0, gamma=7e-6)
    result = fragmentation.simulate(params)
    assert result.population.total_mass == 500
    assert result.time >= 5.0


def test_rings_run():
    params = rings.RingsParams(ntot=300, time_max=1e5)
    result = rings.simulate(params)
    state = result.state
    assert state.chains.total_mass + state.ring_mass == 300
    assert state.n_chains == 1 or result.time >= params.time_max


def test_invalid_parameters():
    with pytest.raises(ValueError):
        fragmentation.FragmentationParams(kd=0.0)
    with pytest.raises(ValueError):
        fragmentation.FragmentationParams(gamma=0.0)
    with pytest.raises(ValueError):
        rings.RingsParams(alpha=1.5)


def test_cli_writes_original_file_names(tmp_path):
    main(["frag", "--ntot", "200", "--time-max", "2", "--gamma", "7e-6", "--outdir", str(tmp_path)])
    names = sorted(p.name for p in tmp_path.iterdir())
    tag = "ntot200_Kd1.00e-03_density1.00e+00_gamma7.00e-06"
    assert f"histograms_{tag}" in names
    assert f"mav_vs_t_{tag}_tmax2.00e+00.dat" in names
    final = [n for n in names if n.startswith(f"masses_final_{tag}_t")]
    assert len(final) == 1
    masses = np.loadtxt(tmp_path / final[0], dtype=int)
    assert masses.shape == (200,) and masses.sum() == 200
    series = np.loadtxt(tmp_path / f"mav_vs_t_{tag}_tmax2.00e+00.dat")
    assert series.shape[1] == 3

    main(["rings", "--ntot", "100", "--outdir", str(tmp_path)])
    assert (tmp_path / "nchain_rings_ntot100_density1.00e+00_pref0.50.dat").exists()
