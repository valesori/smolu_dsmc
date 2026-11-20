"""DSMC of irreversible condensation with ring formation (Ref. 2).

Linear chains anneal with the de Gennes kernel, and a chain of mass ``m``
closes into a ring with rate ``prefactor * m**-1.5``. Rings are inert. The
scheme is the one of :mod:`smolu_dsmc.fragmentation`, with cyclization in
place of fragmentation.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np

from .kernels import cyclization_rate, de_gennes_kernel
from .population import Population, UniformStream


@dataclass(frozen=True)
class RingsParams:
    """Parameters of a simulation; the defaults are those of the original script."""

    ntot: int = 1_000  # Number of monomers
    density: float = 1.0  # Monomer number density
    prefactor: float = 0.5  # Prefactor of the cyclization rate
    alpha: float = 1.0  # Share of the time increment given to annealing, in [0, 1]
    time_max: float = 1.0e5  # Stop when the waiting time exceeds this value
    seed: int | None = 123

    def __post_init__(self) -> None:
        if self.ntot < 2:
            raise ValueError("ntot must be at least 2.")
        if self.density <= 0:
            raise ValueError("density must be > 0.")
        if self.prefactor < 0:
            raise ValueError("prefactor must be >= 0.")
        if not 0.0 <= self.alpha <= 1.0:
            raise ValueError("alpha must be in [0, 1].")


@dataclass
class RingsState:
    """Chains and rings at a given time, with the averages written to output."""

    ntot: int  # Total mass, conserved
    chains: Population
    rings: list[int] = field(default_factory=list)
    ring_mass: int = 0  # Total mass in rings

    @property
    def n_chains(self) -> int:
        return self.chains.n

    @property
    def n_rings(self) -> int:
        return len(self.rings)

    def mean_mass(self) -> float:
        return self.ntot / (self.n_chains + self.n_rings)

    def mean_chain_mass(self) -> float:
        n = self.n_chains
        return ((self.ntot - self.ring_mass) / n) if n else 0.0

    def mean_ring_mass(self) -> float:
        return (self.ring_mass / self.n_rings) if self.rings else 0.0


@dataclass
class RingsResult:
    state: RingsState
    time: float
    steps: int


EventCallback = Callable[[float, RingsState], None]


def simulate(params: RingsParams, *, on_event: EventCallback | None = None) -> RingsResult:
    """Run until the waiting time exceeds ``params.time_max`` or one chain is left."""
    state = RingsState(params.ntot, Population.monomers(params.ntot))
    chains = state.chains
    uniform = UniformStream(np.random.default_rng(params.seed))
    ntot, density, alpha, prefactor = params.ntot, params.density, params.alpha, params.prefactor

    # Exact maxima of both rates, reached by monomers. The annealing one is
    # raised during the run, since the de Gennes kernel grows with mass.
    k_max = de_gennes_kernel(1, 1)
    r_max = cyclization_rate(1, prefactor)

    time = 0.0
    step = 0
    # With a single chain left the annealing probability is 0, and with
    # alpha = 1 the clock would no longer advance: stop there.
    while time < params.time_max and chains.n > 1:
        step += 1
        n = chains.n
        p_anneal = 1.0 / (1.0 + 2.0 * ntot * r_max / ((n - 1) * density * k_max))

        if uniform() < p_anneal:
            i, j = chains.random_pair(uniform)
            k_ij = de_gennes_kernel(chains.mass(i), chains.mass(j))
            if k_ij > k_max:
                k_max = k_ij
            elif uniform() * k_max < k_ij:
                time += 2.0 * alpha * ntot / (n * (n - 1) * density * k_ij)
                chains.merge(i, j)
                if on_event is not None:
                    on_event(time, state)
        else:
            pos = chains.random_molecule(uniform)
            r = cyclization_rate(chains.mass(pos), prefactor)
            if r > r_max:
                r_max = r
            elif uniform() * r_max < r:
                time += (1.0 - alpha) / (n * r)
                mass = chains.remove(pos)
                state.rings.append(mass)
                state.ring_mass += mass
                if on_event is not None:
                    on_event(time, state)

    return RingsResult(state, time, step)
