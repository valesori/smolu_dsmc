"""Bookkeeping of the molecules present in a DSMC simulation.

The original scripts stored masses in a fixed array of length ``ntot`` in which
0 marks an empty slot, and searched it (``np.where``) at every Monte Carlo
step. Here the molecules are kept in a compact list, and the molecules of mass
> 1 in an indexable set, so that every operation used by the algorithm costs
O(1). The sampling probabilities are the same as in the original scripts.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np


class UniformStream:
    """Uniform random numbers in [0, 1), drawn from a generator in blocks.

    Drawing scalars one by one from :class:`numpy.random.Generator` is slow;
    drawing them in blocks gives the same stream several times faster.
    """

    def __init__(self, rng: np.random.Generator, block_size: int = 1 << 16):
        self._rng = rng
        self._block_size = block_size
        self._block: list[float] = []
        self._next = 0

    def __call__(self) -> float:
        if self._next == len(self._block):
            self._block = self._rng.random(self._block_size).tolist()
            self._next = 0
        value = self._block[self._next]
        self._next += 1
        return value

    def index(self, n: int) -> int:
        """Uniform random integer in ``range(n)``."""
        return int(self() * n)


class Population:
    """Masses of the molecules (filaments or chains) present in the system.

    Molecules are addressed by their position ``0 <= p < n``. Positions are not
    stable: removing a molecule moves the last one into its place.
    """

    def __init__(self, masses: Iterable[int]):
        self._masses: list[int] = [int(m) for m in masses if m > 0]
        # Positions of the molecules with mass > 1, and the inverse map.
        self._poly: list[int] = []
        self._poly_index: list[int] = [-1] * len(self._masses)
        for pos, mass in enumerate(self._masses):
            if mass > 1:
                self._poly_add(pos)

    @classmethod
    def monomers(cls, n: int) -> Population:
        return cls([1] * n)

    # --- Queries -------------------------------------------------------------

    @property
    def n(self) -> int:
        """Number of molecules."""
        return len(self._masses)

    @property
    def n_poly(self) -> int:
        """Number of molecules with mass > 1 (the ones that can fragment)."""
        return len(self._poly)

    @property
    def total_mass(self) -> int:
        return sum(self._masses)

    def mass(self, pos: int) -> int:
        return self._masses[pos]

    def masses(self) -> np.ndarray:
        """Masses of all molecules, as a new array."""
        return np.array(self._masses, dtype=np.int64)

    def histogram(self) -> np.ndarray:
        """Two columns: each mass present, and the number of molecules with that mass."""
        values, counts = np.unique(self.masses(), return_counts=True)
        return np.column_stack((values, counts))

    # --- Random selection ----------------------------------------------------

    def random_molecule(self, uniform: UniformStream) -> int:
        return uniform.index(self.n)

    def random_pair(self, uniform: UniformStream) -> tuple[int, int]:
        """A uniformly random ordered pair of distinct molecules."""
        i = uniform.index(self.n)
        j = uniform.index(self.n - 1)
        if j >= i:
            j += 1
        return i, j

    def random_polymer(self, uniform: UniformStream) -> int:
        """A uniformly random molecule among those with mass > 1."""
        return self._poly[uniform.index(self.n_poly)]

    # --- Updates -------------------------------------------------------------

    def merge(self, i: int, j: int) -> None:
        """Replace molecules ``i`` and ``j`` with one molecule of their total mass."""
        self._set_mass(j, self._masses[i] + self._masses[j])
        self.remove(i)

    def split(self, pos: int, m1: int, m2: int) -> None:
        """Break molecule ``pos`` into two molecules of masses ``m1`` and ``m2``."""
        if m1 + m2 != self._masses[pos] or m1 < 1 or m2 < 1:
            raise ValueError(f"cannot split mass {self._masses[pos]} into {m1} + {m2}")
        self._set_mass(pos, m1)
        self.add(m2)

    def add(self, mass: int) -> None:
        self._masses.append(mass)
        self._poly_index.append(-1)
        if mass > 1:
            self._poly_add(self.n - 1)

    def remove(self, pos: int) -> int:
        """Remove molecule ``pos`` and return its mass."""
        mass = self._masses[pos]
        if mass > 1:
            self._poly_remove(pos)
        last = self.n - 1
        if pos != last:
            moved = self._masses[last]
            self._masses[pos] = moved
            k = self._poly_index[last]
            self._poly_index[pos] = k
            if k >= 0:
                self._poly[k] = pos
        self._masses.pop()
        self._poly_index.pop()
        return mass

    def _set_mass(self, pos: int, mass: int) -> None:
        was_poly = self._masses[pos] > 1
        self._masses[pos] = mass
        if mass > 1 and not was_poly:
            self._poly_add(pos)
        elif mass <= 1 and was_poly:
            self._poly_remove(pos)

    def _poly_add(self, pos: int) -> None:
        self._poly_index[pos] = len(self._poly)
        self._poly.append(pos)

    def _poly_remove(self, pos: int) -> None:
        k = self._poly_index[pos]
        last_pos = self._poly[-1]
        self._poly[k] = last_pos
        self._poly_index[last_pos] = k
        self._poly.pop()
        self._poly_index[pos] = -1
