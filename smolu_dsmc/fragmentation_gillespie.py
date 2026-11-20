"""Gillespie simulation of the Smoluchowski equation with fragmentation.

THIS IS A GILLESPIE ALGORITHM, NOT THE DSMC ALGORITHM USED IN REF. 1. The
results of the paper were obtained with :func:`smolu_dsmc.fragmentation.simulate`,
which follows the DSMC scheme of Appendix B of Ref. 1 (after Garcia et al.,
Physica A 143, 535 (1987)). This module samples the same stochastic dynamics
(the same reaction rates, so the same Smoluchowski equation in the
large-system limit) with Gillespie's stochastic simulation algorithm, which is
much faster when filaments are long. It is validated against the DSMC and
against the equilibrium mean length (Ref. 1, Eq. 14) in the tests.

References
----------
- D. T. Gillespie, "Exact stochastic simulation of coupled chemical
  reactions", J. Phys. Chem. 81, 2340-2361 (1977),
  doi:10.1021/j100540a008. The stochastic simulation algorithm.
- D. T. Gillespie, "An exact method for numerically simulating the stochastic
  coalescence process in a cloud", J. Atmos. Sci. 32, 1977-1989 (1975).
  The same method applied to coagulation (stochastic coalescence).
- P. A. W. Lewis and G. S. Shedler, "Simulation of nonhomogeneous Poisson
  processes by thinning", Naval Res. Logist. Q. 26, 403-413 (1979),
  doi:10.1002/nav.3800260304. Thinning: proposing events at an upper-bound
  rate and accepting them with probability rate / bound.
- A. Slepoy, A. P. Thompson and S. J. Plimpton, "A constant-time kinetic Monte
  Carlo algorithm for simulation of large biochemical reaction networks",
  J. Chem. Phys. 128, 205101 (2008), doi:10.1063/1.2919546.
  Composition-rejection selection of weighted events.

Why the algorithm of Appendix B is slow
---------------------------------------
It proposes reactions uniformly and accepts them with probability
``rate / max_rate``. The maximum annealing rate is that of two monomers, while
two filaments of length L anneal ~L times slower; the maximum fragmentation
rate ``(k - 1) F`` grows as ``L_max / 2``, while typical ones are O(1). At
equilibrium more than 95% of the Monte Carlo steps are rejected.

What this module does instead
-----------------------------
It is Gillespie's stochastic simulation algorithm (Gillespie 1975, 1977)
combined with thinning (Lewis and Shedler 1979): reactions are proposed with
probabilities close to their actual rates, accepted with probability
``rate / proposal rate`` (at least ~1/2 here), and the clock advances by an
exponential waiting time drawn from the total proposal rate, rejected
proposals included. It uses the structure of the kernels:

- Annealing. ``K_ij = (1/i + 1/j) / K_d`` is a sum of a term depending on
  ``i`` only and one depending on ``j`` only. Picking filament ``a`` with
  probability proportional to ``1/m_a`` and ``b`` uniformly among the others
  gives each pair a probability proportional to ``K_ab``. The total rate is
  ``(n - 1) sum_a 1/m_a / (K_d V)``.
- Fragmentation. The rate at which a k-mer breaks into ``(j, k - j)`` is
  ``(1/j + 1/(k - j)) / 2`` in the units of Ref. 1 (Eq. 1 and 16), so its
  total fragmentation rate is the harmonic number ``H_{k-1} / 2``. Picking the
  filament with probability proportional to ``H_{k-1}``, and then ``j`` with
  probability proportional to ``1/j``, samples each fragmentation with the
  right probability.
- Weighted choices use composition-rejection (Slepoy et al. 2008): filaments
  are grouped by
  ``floor(log2(mass))``, a group is chosen with probability proportional to
  an upper bound of its total weight, and a member is chosen uniformly and
  accepted with probability ``weight / bound`` (>= 1/2).
- Entanglement. The free kernel is used for proposals, and a proposal is
  accepted with an additional probability ``g(x, gamma) <= 1``.

Differences in behaviour from the DSMC of :func:`smolu_dsmc.fragmentation.simulate`:
waiting times are exponential random variables rather than their mean value
(so there is no ``alpha`` parameter, and single trajectories, not only
averages, follow the exact stochastic dynamics), and a "step" (used for
histogram snapshots) is one proposal, of which a much larger fraction is
accepted.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence

import numpy as np

from .fragmentation import FragmentationParams, FragmentationResult
from .kernels import blocking_factor
from .population import UniformStream
from .recording import RecordCallback, TimeGrid

_EULER_GAMMA = 0.5772156649015329
_RECOMPUTE_EVERY = 1 << 16  # Events between full recomputations of the bounds


def _harmonic_table(size: int) -> list[float]:
    table = [0.0] * size
    for i in range(1, size):
        table[i] = table[i - 1] + 1.0 / i
    return table


_HARMONIC_TABLE = _harmonic_table(1 << 12)


def harmonic(n: int) -> float:
    """Harmonic number ``H_n = 1 + 1/2 + ... + 1/n`` (``H_0 = 0``)."""
    if n < len(_HARMONIC_TABLE):
        return _HARMONIC_TABLE[n]
    inv = 1.0 / n
    inv2 = inv * inv
    # Asymptotic expansion; relative error below 1e-16 for n >= 4096.
    return math.log(n) + _EULER_GAMMA + 0.5 * inv - inv2 / 12 + inv2 * inv2 / 120


def sample_inverse_index(k: int, uniform: UniformStream) -> int:
    """Random ``j`` in ``1..k`` with probability proportional to ``1/j``.

    Proposes ``j = floor(x)`` with ``x`` from the density ``1/x`` on
    ``[1, k+1)``, which gives ``P(j)`` proportional to ``ln(1 + 1/j)``, and
    accepts with probability ``ln(2) / (j ln(1 + 1/j))`` (between ~0.69 and 1).
    """
    log_top = math.log(k + 1)
    while True:
        j = int(math.exp(uniform() * log_top))
        if j > k:  # Possible only through rounding.
            continue
        if uniform() * j * math.log1p(1.0 / j) < math.log(2.0):
            return j


class GroupedPopulation:
    """Molecules grouped by ``floor(log2(mass))`` for weighted random choices.

    Offers the same read-only interface as
    :class:`smolu_dsmc.population.Population` (``n``, ``total_mass``,
    ``mass``, ``masses``, ``histogram``), which the output writers use.
    """

    def __init__(self, masses, max_mass: int):
        n_groups = max(1, int(max_mass).bit_length())
        self._masses: list[int] = []
        self._slot: list[int] = []  # Index of each molecule in its group list
        self.groups: list[list[int]] = [[] for _ in range(n_groups)]
        # Upper bounds, for members of group g, of 1/m and of H_{m-1}.
        self.inv_bound = [2.0**-g for g in range(n_groups)]
        self.harm_bound = [harmonic(2 ** (g + 1) - 2) for g in range(n_groups)]
        self.inv_total = 0.0
        self.harm_total = 0.0
        for m in masses:
            if m > 0:
                self.add(int(m))
        self.recompute_totals()

    # --- Read-only interface shared with Population --------------------------

    @property
    def n(self) -> int:
        return len(self._masses)

    @property
    def total_mass(self) -> int:
        return sum(self._masses)

    def mass(self, pos: int) -> int:
        return self._masses[pos]

    def masses(self) -> np.ndarray:
        return np.array(self._masses, dtype=np.int64)

    def histogram(self) -> np.ndarray:
        values, counts = np.unique(self.masses(), return_counts=True)
        return np.column_stack((values, counts))

    # --- Bounds on the total weights -----------------------------------------

    def recompute_totals(self) -> None:
        """Sums over molecules of the bounds on ``1/m`` and on ``H_{m-1}``."""
        self.inv_total = math.fsum(
            len(g) * w for g, w in zip(self.groups, self.inv_bound, strict=True)
        )
        self.harm_total = math.fsum(
            len(g) * w for g, w in zip(self.groups, self.harm_bound, strict=True)
        )

    def _choose(self, weights: list[float], total: float, uniform: UniformStream) -> int:
        """A molecule chosen with probability proportional to its group bound."""
        target = uniform() * total
        groups = self.groups
        last = 0
        for g, members in enumerate(groups):
            if members:
                last = g
                target -= len(members) * weights[g]
                if target < 0.0:
                    return members[uniform.index(len(members))]
        members = groups[last]  # Rounding: fall back on the last non-empty group.
        return members[uniform.index(len(members))]

    def choose_by_inverse_mass(self, uniform: UniformStream) -> int | None:
        """A molecule with probability proportional to ``1/m``, or None (rejected)."""
        pos = self._choose(self.inv_bound, self.inv_total, uniform)
        m = self._masses[pos]
        if uniform() * self.inv_bound[m.bit_length() - 1] * m < 1.0:
            return pos
        return None

    def choose_by_harmonic(self, uniform: UniformStream) -> int | None:
        """A molecule with probability proportional to ``H_{m-1}``, or None (rejected)."""
        pos = self._choose(self.harm_bound, self.harm_total, uniform)
        m = self._masses[pos]
        if uniform() * self.harm_bound[m.bit_length() - 1] < harmonic(m - 1):
            return pos
        return None

    # --- Updates -------------------------------------------------------------

    def add(self, mass: int) -> None:
        pos = len(self._masses)
        self._masses.append(mass)
        self._slot.append(0)
        self._group_add(pos, mass)

    def remove(self, pos: int) -> int:
        mass = self._masses[pos]
        self._group_remove(pos, mass)
        last = len(self._masses) - 1
        if pos != last:
            moved = self._masses[last]
            g = self.groups[moved.bit_length() - 1]
            slot = self._slot[last]
            g[slot] = pos
            self._masses[pos] = moved
            self._slot[pos] = slot
        self._masses.pop()
        self._slot.pop()
        return mass

    def set_mass(self, pos: int, mass: int) -> None:
        old = self._masses[pos]
        if old.bit_length() != mass.bit_length():
            self._group_remove(pos, old)
            self._group_add(pos, mass)
        self._masses[pos] = mass

    def _group_add(self, pos: int, mass: int) -> None:
        g = mass.bit_length() - 1
        members = self.groups[g]
        self._slot[pos] = len(members)
        members.append(pos)
        self.inv_total += self.inv_bound[g]
        self.harm_total += self.harm_bound[g]

    def _group_remove(self, pos: int, mass: int) -> None:
        g = mass.bit_length() - 1
        members = self.groups[g]
        slot = self._slot[pos]
        last_pos = members[-1]
        members[slot] = last_pos
        self._slot[last_pos] = slot
        members.pop()
        self.inv_total -= self.inv_bound[g]
        self.harm_total -= self.harm_bound[g]


EventCallback = Callable[[float, GroupedPopulation], None]
SnapshotCallback = Callable[[int, float, GroupedPopulation], None]


def simulate_gillespie(
    params: FragmentationParams,
    *,
    masses=None,
    on_event: EventCallback | None = None,
    on_snapshot: SnapshotCallback | None = None,
    snapshot_every: int = 10_000,
    record_times: Sequence[float] | None = None,
    on_record: RecordCallback | None = None,
) -> FragmentationResult:
    """Gillespie alternative to the DSMC of :func:`smolu_dsmc.fragmentation.simulate`.

    NOT the algorithm used in Ref. 1 (see the module docstring for the method
    and references); it samples the same dynamics, faster. ``params.alpha``
    is ignored. ``masses`` defaults to ``params.ntot`` monomers. The callbacks
    are those of :func:`smolu_dsmc.fragmentation.simulate`; here a "step" is
    one proposal.
    """
    if masses is None:
        masses = [1] * params.ntot
    pop = GroupedPopulation(masses, max_mass=params.ntot)
    ntot = pop.total_mass
    if ntot != params.ntot:
        raise ValueError(f"masses have total {ntot}, expected {params.ntot}")

    uniform = UniformStream(np.random.default_rng(params.seed))
    density, kd, gamma = params.density, params.kd, params.gamma
    # Annealing rate per unordered pair is K_ab / V, with V = ntot / density.
    anneal_prefactor = density / (ntot * kd)

    grid = TimeGrid(record_times, on_record) if on_record and record_times is not None else None
    next_record = grid.next_time if grid else math.inf

    time = 0.0
    step = 0
    events = 0
    while time < params.time_max:
        if on_snapshot is not None and step % snapshot_every == 0:
            on_snapshot(step, time, pop)
        step += 1

        n = pop.n
        anneal_bound = anneal_prefactor * (n - 1) * pop.inv_total
        frag_bound = 0.5 * pop.harm_total
        total = anneal_bound + frag_bound
        if total <= 0.0:
            break  # A single monomer: nothing can happen.
        time += -math.log(1.0 - uniform()) / total
        if time > next_record:
            next_record = grid.record_until(time, pop)
        if time >= params.time_max:
            break

        if uniform() * total < anneal_bound:
            a = pop.choose_by_inverse_mass(uniform)
            if a is None:
                continue
            b = uniform.index(n - 1)
            if b >= a:
                b += 1
            ma, mb = pop.mass(a), pop.mass(b)
            if gamma is not None:
                l_min = min(ma, mb)
                if uniform() >= blocking_factor(0.5 * density * l_min * l_min, gamma):
                    continue
            pop.set_mass(b, ma + mb)
            pop.remove(a)
        else:
            pos = pop.choose_by_harmonic(uniform)
            if pos is None:
                continue
            mk = pop.mass(pos)
            m1 = sample_inverse_index(mk - 1, uniform)
            m2 = mk - m1
            if gamma is not None:
                l_min = min(m1, m2)
                if uniform() >= blocking_factor(0.5 * density * l_min * l_min, gamma):
                    continue
            pop.set_mass(pos, m1)
            pop.add(m2)

        events += 1
        if events % _RECOMPUTE_EVERY == 0:
            pop.recompute_totals()
        if on_event is not None:
            on_event(time, pop)

    if grid is not None:
        grid.record_remaining(pop)
    return FragmentationResult(pop, time, step)
