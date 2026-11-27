"""DSMC of the Smoluchowski equation with annealing and fragmentation (Ref. 1).

The algorithm is the one described in Appendix B of Ref. 1, which follows
Garcia et al., Physica A 143, 535 (1987). Each Monte Carlo step attempts an
annealing with probability ``p`` (Eq. 17) or a fragmentation with probability
``1 - p``; the attempt is accepted with probability ``rate / rate_max``. The
waiting time advances only on accepted events, by ``alpha`` times the annealing
increment and ``1 - alpha`` times the fragmentation increment, which gives the
correct kinetics on average for any ``0 <= alpha <= 1``.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from functools import partial
from typing import Any

import numpy as np

from .kernels import entangled_kernel, rouse_kernel
from .population import Population, UniformStream
from .recording import RecordCallback, TimeGrid

Kernel = Callable[[int, int], float]

AVOGADRO = 6.02214076e23  # 1/mol
# Reference density c0 of Ref. 1: the lowest experimental vimentin concentration,
# 0.01 mg/mL, i.e. 6.1 nM of ULFs (K_d = 1e-3 c0 = 1.0e-5 mg/mL = 6.1 pM, Ref. 1,
# Sec. II.B), in ULFs per m^3 (about 3.67e18).
REF1_REFERENCE_DENSITY = 6.1e-9 * 1e3 * AVOGADRO
# Axial length per ULF in Ref. 1 (nominal ULF length), in m.
REF1_ULF_AXIAL_LENGTH = 60e-9

# Values of FragmentationParams.entanglement_x (see its docstring).
ENTANGLEMENT_X = ("physical", "reduced")


@dataclass(frozen=True)
class FragmentationParams:
    """Parameters of a simulation; the defaults are those of the original scripts.

    With entanglement (``gamma`` not None), ``entanglement_x`` selects how the
    variable ``x`` of the blocking factor is computed (see the README,
    "Entanglement"):

    - ``"physical"`` (default): from the physical density,
      ``x = c0 a**3 * density * l_min**2 / 2``, and ``gamma`` is the physical
      maximum misalignment angle, in rad. ``reference_density`` (c0, in m^-3)
      and ``axial_length`` (a, axial length per ULF, in m) default to the
      values of Ref. 1.
    - ``"reduced"``: the legacy convention, from the reduced density,
      ``x = density * l_min**2 / 2``, as in the 2022 code (``legacy_2022/``);
      ``gamma`` is then the angle fitted in Ref. 1 (7e-6).
      ``reference_density`` and ``axial_length`` must be left unset.

    In both, ``density`` is the reduced density c/c0, used unchanged by the
    rest of the simulation.
    """

    ntot: int = 10_000  # Number of monomers (ULFs)
    kd: float = 1.0e-3  # Equilibrium dissociation constant K_d, in units of c0
    density: float = 1.0  # Reduced ULF number density c/c0
    alpha: float = 1.0  # Share of the time increment given to annealing, in [0, 1]
    time_max: float = 1.0e5  # Stop when the waiting time exceeds this value
    gamma: float | None = None  # Max. annealing angle in rad; None disables entanglement
    seed: int | None = 123
    entanglement_x: str = "physical"  # "physical" or "reduced" (legacy)
    # entanglement_x == "physical" only; None selects the Ref. 1 values above.
    reference_density: float | None = None  # c0, in m^-3
    axial_length: float | None = None  # a, axial length per ULF, in m

    def __post_init__(self) -> None:
        if self.ntot < 2:
            raise ValueError("ntot must be at least 2.")
        if self.kd <= 0:
            raise ValueError("kd must be > 0.")
        if self.density <= 0:
            raise ValueError("density must be > 0.")
        if not 0.0 <= self.alpha <= 1.0:
            raise ValueError("alpha must be in [0, 1].")
        if self.gamma is not None and not 0.0 < self.gamma <= math.pi / 2:
            raise ValueError("gamma must be in (0, pi/2].")
        if self.entanglement_x not in ENTANGLEMENT_X:
            raise ValueError(
                f"entanglement_x must be one of {ENTANGLEMENT_X}, not {self.entanglement_x!r}."
            )
        physical = {"reference_density": self.reference_density, "axial_length": self.axial_length}
        if self.gamma is not None and self.entanglement_x == "physical":
            defaults = {
                "reference_density": REF1_REFERENCE_DENSITY,
                "axial_length": REF1_ULF_AXIAL_LENGTH,
            }
            for name, value in physical.items():
                if value is None:
                    object.__setattr__(self, name, defaults[name])
                elif not value > 0:
                    raise ValueError(f"{name} must be > 0.")
        else:
            given = [name for name, value in physical.items() if value is not None]
            if given:
                raise ValueError(
                    f"{' and '.join(given)} only apply with gamma set and "
                    "entanglement_x = 'physical'."
                )

    @property
    def legacy_entanglement(self) -> bool:
        """True for an entangled run with the legacy (reduced-density) x."""
        return self.gamma is not None and self.entanglement_x == "reduced"

    @property
    def entanglement_c0_a3(self) -> float:
        """Factor converting ``density`` into the density in ``x``: c0 a^3, or 1 (reduced)."""
        if self.gamma is None or self.entanglement_x == "reduced":
            return 1.0
        return self.reference_density * self.axial_length**3

    def entanglement_record(self) -> dict[str, Any] | None:
        """Description of the entanglement variable, for the run record."""
        if self.gamma is None:
            return None
        if self.entanglement_x == "reduced":
            return {
                "x": "reduced",
                "x_formula": "density * l_min**2 / 2",
                "gamma_rad": self.gamma,
                "gamma_meaning": "fitted angle of the legacy convention (Ref. 1 fit: 7e-6)",
            }
        return {
            "x": "physical",
            "x_formula": "c0_a3 * density * l_min**2 / 2",
            "gamma_rad": self.gamma,
            "gamma_meaning": "physical maximum misalignment angle",
            "reference_density_per_m3": self.reference_density,
            "axial_length_m": self.axial_length,
            "c0_a3": self.entanglement_c0_a3,
        }

    @property
    def kernel(self) -> Kernel:
        if self.gamma is None:
            return rouse_kernel
        return partial(
            entangled_kernel,
            density=self.density,
            gamma=self.gamma,
            c0_a3=self.entanglement_c0_a3,
        )


@dataclass
class FragmentationResult:
    population: Population
    time: float
    steps: int
    # Final estimates of the maximum rates (Appendix B algorithm only).
    k_max: float | None = None
    f_max: float | None = None


# Called after each accepted event with (time, population).
EventCallback = Callable[[float, Population], None]
# Called every ``snapshot_every`` Monte Carlo steps with (step, time, population).
SnapshotCallback = Callable[[int, float, Population], None]


def simulate(
    params: FragmentationParams,
    *,
    population: Population | None = None,
    on_event: EventCallback | None = None,
    on_snapshot: SnapshotCallback | None = None,
    snapshot_every: int = 10_000,
    record_times: Sequence[float] | None = None,
    on_record: RecordCallback | None = None,
) -> FragmentationResult:
    """Run the simulation until the waiting time exceeds ``params.time_max``.

    ``population`` defaults to ``params.ntot`` monomers. Three optional ways
    to observe the run:

    - ``on_event(time, population)`` after each accepted event;
    - ``on_snapshot(step, time, population)`` every ``snapshot_every`` Monte
      Carlo steps;
    - ``on_record(index, t, population)`` with the state at each time ``t`` of
      ``record_times`` (see :mod:`smolu_dsmc.recording`).
    """
    pop = population if population is not None else Population.monomers(params.ntot)
    ntot = pop.total_mass
    if ntot != params.ntot:
        raise ValueError(f"population has total mass {ntot}, expected {params.ntot}")

    uniform = UniformStream(np.random.default_rng(params.seed))
    kernel = params.kernel
    kd, density, alpha = params.kd, params.density, params.alpha

    # Initial estimates of the maximum rates; both are raised when exceeded.
    # k_max is the exact maximum of K_ij; f_max that of (k - 1) F_{i,k-i} for k = 2.
    k_max = kernel(1, 1) / kd
    f_max = 0.5 * kernel(1, 1)

    grid = TimeGrid(record_times, on_record) if on_record and record_times is not None else None
    next_record = grid.next_time if grid else math.inf

    time = 0.0
    step = 0
    while time < params.time_max:
        if on_snapshot is not None and step % snapshot_every == 0:
            on_snapshot(step, time, pop)
        step += 1

        n, n_poly = pop.n, pop.n_poly
        if n < 2:
            if n_poly == 0:
                break  # A single monomer: nothing can happen.
            p_anneal = 0.0
        else:
            p_anneal = 1.0 / (1.0 + ntot * n_poly * f_max / (n * (n - 1) * density * k_max))

        if uniform() < p_anneal:
            i, j = pop.random_pair(uniform)
            k_ij = kernel(pop.mass(i), pop.mass(j)) / kd
            if k_ij > k_max:
                k_max = k_ij  # Raise the estimate and discard the attempt (Ref. 1, App. B).
            elif uniform() * k_max < k_ij:
                time += 2.0 * alpha * ntot / (n * (n - 1) * density * k_ij)
                if time > next_record:
                    next_record = grid.record_until(time, pop)
                pop.merge(i, j)
                if on_event is not None:
                    on_event(time, pop)
        elif n_poly > 0:
            pos = pop.random_polymer(uniform)
            mk = pop.mass(pos)
            m1 = 1 + uniform.index(mk - 1)
            m2 = mk - m1
            # Total rate (k - 1) F, so that one split chosen uniformly is accepted
            # with probability proportional to F.
            f_tot = (mk - 1) * 0.5 * kernel(m1, m2)
            if f_tot > f_max:
                f_max = f_tot
            elif uniform() * f_max < f_tot:
                time += 2.0 * (1.0 - alpha) / (n_poly * f_tot)
                if time > next_record:
                    next_record = grid.record_until(time, pop)
                pop.split(pos, m1, m2)
                if on_event is not None:
                    on_event(time, pop)

    if grid is not None:
        grid.record_remaining(pop)
    return FragmentationResult(pop, time, step, k_max, f_max)
