"""Recording the state of a simulation at prescribed times.

The simulations change state only at discrete events. The state "at time t"
is the one holding between the last event before ``t`` and the first event
after it, so a :class:`TimeGrid` records a grid time ``t`` just before applying
the first event whose time exceeds ``t``.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from typing import Any

import numpy as np

# Called with (index of the grid time, grid time, population).
RecordCallback = Callable[[int, float, Any], None]


def log_time_grid(t_min: float, t_max: float, n_points: int) -> np.ndarray:
    """``n_points`` logarithmically spaced times from ``t_min`` to ``t_max``."""
    if n_points < 2:
        raise ValueError("the time grid needs at least 2 points.")
    if not 0.0 < t_min < t_max:
        raise ValueError(f"the time grid needs 0 < t_min < t_max (got {t_min}, {t_max}).")
    return np.geomspace(t_min, t_max, n_points)


class TimeGrid:
    """Calls ``callback(index, t, population)`` for each grid time ``t``.

    The simulation calls :meth:`record_until` with the time of each event
    before applying it, and :meth:`record_remaining` when it stops.
    """

    def __init__(self, times: Sequence[float], callback: RecordCallback):
        self._times = [float(t) for t in times]
        if any(t < 0 for t in self._times) or any(
            a > b for a, b in zip(self._times, self._times[1:], strict=False)
        ):
            raise ValueError("record times must be non-negative and sorted.")
        self._callback = callback
        self._next = 0

    @property
    def next_time(self) -> float:
        """The next grid time to record, or infinity if none is left."""
        return self._times[self._next] if self._next < len(self._times) else math.inf

    def record_until(self, time: float, population: Any) -> float:
        """Record the grid times before ``time``; return :attr:`next_time`."""
        times = self._times
        while self._next < len(times) and times[self._next] < time:
            self._callback(self._next, times[self._next], population)
            self._next += 1
        return self.next_time

    def record_remaining(self, population: Any) -> None:
        """Record the grid times left when the simulation stops.

        The state no longer changes at that point: either the end time has
        been reached, or no reaction is possible.
        """
        while self._next < len(self._times):
            self._callback(self._next, self._times[self._next], population)
            self._next += 1
