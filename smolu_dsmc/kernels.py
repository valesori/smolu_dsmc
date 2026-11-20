"""Reaction kernels, in the dimensionless units of Refs. 1 and 2.

Annealing and fragmentation (Ref. 1, Eqs. 2-4 and 16):

    K_ij = kernel(i, j) / K_d       (annealing)
    F_ij = kernel(i, j) / 2         (fragmentation, detailed balance F_ij = K_d K_ij / 2)

where ``kernel`` is :func:`rouse_kernel` for freely diffusing filaments, or
:func:`entangled_kernel` when entanglements block part of the reactions.

Ring formation (Ref. 2): :func:`de_gennes_kernel` for annealing and
:func:`cyclization_rate` for loop closure.
"""

from __future__ import annotations

import math


def rouse_kernel(i: int, j: int) -> float:
    """Diffusion-limited end-to-end annealing of Rouse filaments, ``1/i + 1/j``."""
    return 1.0 / i + 1.0 / j


def blocking_factor(x: float, gamma: float) -> float:
    """Fraction of annealing reactions not blocked by entanglements, ``g(x, gamma)``.

    ``g(x, gamma) = [1 - (cos(gamma) + x sin(gamma)) exp(-gamma x)]
                    / [(1 + x^2) (1 - cos(gamma))]``

    (Ref. 1, Sec. II.C). It goes from 1 at ``x = 0`` to 0 as ``x`` grows.

    For the small angles of interest (gamma ~ 1e-5) numerator and denominator
    are ~ gamma**2, and evaluating them as written loses about half of the
    significant digits to cancellation. They are rewritten here as sums of
    positive terms:

        1 - cos(g)                        = 2 sin(g/2)**2
        1 - (cos(g) + x sin(g)) exp(-u)   = exp(-u) [(e^u - 1 - u) + (1 - cos(g)) + x (g - sin(g))]

    with ``u = gamma x``. The second form is used when ``u`` is small, where the
    cancellation happens.
    """
    u = gamma * x
    one_minus_cos = 2.0 * math.sin(0.5 * gamma) ** 2
    if u < 1e-2:
        exp_rest = u * u * (1 / 2 + u * (1 / 6 + u * (1 / 24 + u * (1 / 120 + u / 720))))
        numerator = math.exp(-u) * (exp_rest + one_minus_cos + x * _gamma_minus_sin(gamma))
    else:
        numerator = 1.0 - (math.cos(gamma) + x * math.sin(gamma)) * math.exp(-u)
    return numerator / ((1.0 + x * x) * one_minus_cos)


def _gamma_minus_sin(g: float) -> float:
    """``g - sin(g)``, accurate also for small ``g``."""
    if g < 1e-2:
        g2 = g * g
        return g * g2 / 6 * (1 - g2 / 20 * (1 - g2 / 42 * (1 - g2 / 72)))
    return g - math.sin(g)


def entangled_kernel(i: int, j: int, density: float, gamma: float) -> float:
    """:func:`rouse_kernel` multiplied by the entanglement blocking factor.

    ``x = density * l_min**2 / 2``, with ``l_min = min(i, j)`` (Ref. 1, Eq. 4).
    """
    l_min = min(i, j)
    x = 0.5 * density * l_min * l_min
    return blocking_factor(x, gamma) * rouse_kernel(i, j)


def de_gennes_kernel(i: int, j: int) -> float:
    """Annealing kernel for Gaussian chains, ``(sqrt(i) + sqrt(j)) (1/i + 1/j)``."""
    return (math.sqrt(i) + math.sqrt(j)) * (1.0 / i + 1.0 / j)


def cyclization_rate(m: int, prefactor: float) -> float:
    """Loop-closure rate of a chain of mass ``m``, ``prefactor * m**-1.5``."""
    return prefactor * m**-1.5
