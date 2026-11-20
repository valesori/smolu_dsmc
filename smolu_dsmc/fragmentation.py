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
from collections.abc import Callable
from dataclasses import dataclass
from functools import partial

import numpy as np

from .kernels import entangled_kernel, rouse_kernel
from .population import Population, UniformStream

Kernel = Callable[[int, int], float]


@dataclass(frozen=True)
class FragmentationParams:
    """Parameters of a simulation; the defaults are those of the original scripts."""

    ntot: int = 10_000  # Number of monomers (ULFs)
    kd: float = 1.0e-3  # Equilibrium dissociation constant K_d
    density: float = 1.0  # Monomer number density c
    alpha: float = 1.0  # Share of the time increment given to annealing, in [0, 1]
    time_max: float = 1.0e5  # Stop when the waiting time exceeds this value
    gamma: float | None = None  # Max. annealing angle; None disables entanglement
    seed: int | None = 123

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

    @property
    def kernel(self) -> Kernel:
        if self.gamma is None:
            return rouse_kernel
        return partial(entangled_kernel, density=self.density, gamma=self.gamma)


@dataclass
class FragmentationResult:
    population: Population
    time: float
    steps: int
    k_max: float
    f_max: float


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
) -> FragmentationResult:
    """Run the simulation until the waiting time exceeds ``params.time_max``.

    ``population`` defaults to ``params.ntot`` monomers.
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
                pop.split(pos, m1, m2)
                if on_event is not None:
                    on_event(time, pop)

    return FragmentationResult(pop, time, step, k_max, f_max)
