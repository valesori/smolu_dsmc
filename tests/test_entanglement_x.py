"""Physical (default) and reduced (legacy) entanglement variable x."""

import json
import math

import pytest

from smolu_dsmc import fragmentation
from smolu_dsmc.cli import main
from smolu_dsmc.fragmentation import FragmentationParams, simulate
from smolu_dsmc.fragmentation_gillespie import simulate_gillespie
from smolu_dsmc.kernels import blocking_factor, entanglement_variable, rouse_kernel

AVOGADRO = 6.02214076e23
C0 = 6.1e-9 * 1e3 * AVOGADRO  # 6.1 nM in m^-3, computed independently of the package
A = 60e-9  # m
GAMMA_FIT = 7e-6  # Ref. 1, legacy convention (reduced x)


def test_reference_multiplier():
    params = FragmentationParams(gamma=1e-2)
    assert params.reference_density == pytest.approx(3.67e18, rel=1e-3)
    assert params.axial_length == 60e-9
    assert params.entanglement_c0_a3 == pytest.approx(C0 * A**3, rel=1e-12)
    assert params.entanglement_c0_a3 == pytest.approx(7.9e-4, rel=5e-3)


@pytest.mark.parametrize("density", [1.0, 2.5])
@pytest.mark.parametrize("l_min", [1, 37, 1000])
def test_physical_x(density, l_min):
    """x = c_phys a^3 l_min^2 / 2 with c_phys = density * c0, in SI units."""
    params = FragmentationParams(density=density, gamma=1e-2)
    x = entanglement_variable(l_min, density, params.entanglement_c0_a3)
    c_phys = density * C0  # m^-3
    assert x == pytest.approx(c_phys * A**3 * l_min**2 / 2, rel=1e-12)
    # and x_physical = c0 a^3 x_reduced
    assert x == pytest.approx(C0 * A**3 * density * l_min**2 / 2, rel=1e-12)


@pytest.mark.parametrize("density", [1.0, 2.5])
def test_physical_kernel(density):
    params = FragmentationParams(density=density, gamma=1e-2)
    for i, j in [(1, 1), (40, 300), (5000, 70)]:
        x = density * C0 * A**3 * min(i, j) ** 2 / 2
        expected = blocking_factor(x, 1e-2) * rouse_kernel(i, j)
        assert params.kernel(i, j) == pytest.approx(expected, rel=1e-12)


def test_custom_reference_density_and_axial_length():
    params = FragmentationParams(gamma=1e-2, reference_density=2e18, axial_length=49e-9)
    assert params.entanglement_c0_a3 == pytest.approx(2e18 * (49e-9) ** 3, rel=1e-12)


@pytest.mark.parametrize("density", [1.0, 5.0])
def test_reduced_kernel_is_unchanged(density):
    params = FragmentationParams(density=density, gamma=GAMMA_FIT, entanglement_x="reduced")
    assert params.entanglement_c0_a3 == 1.0
    assert params.reference_density is None and params.axial_length is None
    for i, j in [(1, 1), (40, 300), (5000, 70)]:
        l_min = min(i, j)
        # The expression of the code before the correction, bit for bit.
        expected = blocking_factor(0.5 * density * l_min * l_min, GAMMA_FIT) * rouse_kernel(i, j)
        assert params.kernel(i, j) == expected


def test_fitted_angle_and_physical_angle_block_alike():
    """Small-angle equivalence: g(x_reduced, gamma_fit) ~ g(x_physical, gamma_fit / c0 a^3).

    The difference is at most about 5e-7 (at l_min ~ 700, gamma x ~ 1.8), so
    the tolerance is 1e-6 in absolute value.
    """
    c0_a3 = C0 * A**3
    gamma_phys = GAMMA_FIT / c0_a3
    assert gamma_phys == pytest.approx(8.82e-3, rel=1e-3)
    for density in (1.0, 5.0):
        for l_min in [*range(1, 3000), *range(3000, 100_001, 97)]:
            g_legacy = blocking_factor(entanglement_variable(l_min, density, 1.0), GAMMA_FIT)
            g_physical = blocking_factor(entanglement_variable(l_min, density, c0_a3), gamma_phys)
            assert g_physical == pytest.approx(g_legacy, abs=1e-6)


def test_same_angle_differs_between_x():
    """The angle is not rescaled: the same gamma means different blocking."""
    legacy = FragmentationParams(gamma=GAMMA_FIT, entanglement_x="reduced")
    physical = FragmentationParams(gamma=GAMMA_FIT)
    assert physical.gamma == legacy.gamma == GAMMA_FIT
    assert legacy.kernel(1000, 1000) < 0.5 * physical.kernel(1000, 1000)


@pytest.mark.parametrize("simulate_fn", [simulate, simulate_gillespie])
def test_x_irrelevant_without_entanglement(simulate_fn):
    results = [
        simulate_fn(FragmentationParams(ntot=500, kd=1e-2, time_max=5.0, seed=3, entanglement_x=x))
        for x in fragmentation.ENTANGLEMENT_X
    ]
    a, b = results
    assert (a.time, a.steps) == (b.time, b.steps)
    assert sorted(a.population.masses()) == sorted(b.population.masses())


@pytest.mark.parametrize("simulate_fn", [simulate, simulate_gillespie])
def test_physical_run_conserves_mass(simulate_fn):
    result = simulate_fn(FragmentationParams(ntot=500, time_max=5.0, gamma=8.8e-3))
    assert result.population.total_mass == 500
    assert result.time >= 5.0


def test_invalid_entanglement_parameters():
    with pytest.raises(ValueError, match="entanglement_x"):
        FragmentationParams(gamma=1e-2, entanglement_x="legacy")
    with pytest.raises(ValueError, match="only apply"):
        FragmentationParams(gamma=7e-6, entanglement_x="reduced", axial_length=60e-9)
    with pytest.raises(ValueError, match="only apply"):
        FragmentationParams(reference_density=3e18)  # no gamma
    with pytest.raises(ValueError, match="axial_length"):
        FragmentationParams(gamma=1e-2, axial_length=-1.0)


def test_cli_physical_run_record(tmp_path):
    main(
        [
            "frag",
            "--ntot",
            "200",
            "--time-max",
            "2",
            "--gamma",
            "8.82e-3",
            "--outdir",
            str(tmp_path),
        ]
    )
    tag = "ntot200_Kd1.00e-03_density1.00e+00_gamma8.82e-03"
    assert (tmp_path / f"mav_vs_t_{tag}_tmax2.00e+00.dat").exists()
    record = json.loads((tmp_path / f"run_{tag}_tmax2.00e+00.json").read_text())
    ent = record["entanglement"]
    assert ent["x"] == "physical"
    assert ent["gamma_rad"] == 8.82e-3
    assert ent["reference_density_per_m3"] == pytest.approx(C0)
    assert ent["axial_length_m"] == 60e-9
    assert ent["c0_a3"] == pytest.approx(C0 * A**3)
    assert record["parameters"]["entanglement_x"] == "physical"


def test_cli_legacy_run_record(tmp_path):
    main(
        [
            "frag",
            "--ntot",
            "200",
            "--time-max",
            "2",
            "--gamma",
            "7e-6",
            "--entanglement-x",
            "reduced",
            "--outdir",
            str(tmp_path),
        ]
    )
    tag = "ntot200_Kd1.00e-03_density1.00e+00_gamma7.00e-06_legacy"
    record = json.loads((tmp_path / f"run_{tag}_tmax2.00e+00.json").read_text())
    assert record["entanglement"]["x"] == "reduced"
    assert record["entanglement"]["gamma_rad"] == 7e-6
    assert "reference_density_per_m3" not in record["entanglement"]
    assert record["parameters"]["reference_density"] is None


def test_cli_free_run_record(tmp_path):
    main(["frag", "--ntot", "200", "--time-max", "2", "--outdir", str(tmp_path)])
    record = json.loads(
        (tmp_path / "run_ntot200_Kd1.00e-03_density1.00e+00_tmax2.00e+00.json").read_text()
    )
    assert record["entanglement"] is None


def test_cli_rejects_physical_constants_with_reduced_x(tmp_path):
    with pytest.raises(SystemExit):
        main(
            [
                "frag",
                "--gamma",
                "7e-6",
                "--entanglement-x",
                "reduced",
                "--axial-length",
                "49e-9",
                "--outdir",
                str(tmp_path),
            ]
        )


def test_physical_angle_of_ref1_fit():
    """The value quoted in the README and the CLI help."""
    assert math.isclose(GAMMA_FIT / (C0 * A**3), 8.822e-3, rel_tol=1e-3)
