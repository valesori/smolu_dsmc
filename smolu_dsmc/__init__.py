"""Stochastic simulation of Smoluchowski equations.

Models:

- :mod:`smolu_dsmc.fragmentation`: annealing and fragmentation, optionally with
  entanglement, simulated with the Direct Simulation Monte Carlo (DSMC)
  algorithm of Ref. 1 (Appendix B), which was used for its published results.
- :mod:`smolu_dsmc.fragmentation_gillespie`: the same model simulated with a
  Gillespie algorithm. It samples the same dynamics faster, but it is not the
  algorithm used in Ref. 1.
- :mod:`smolu_dsmc.rings`: irreversible condensation with ring formation
  (Ref. 2), simulated with DSMC.

Quantities are dimensionless; see the README for the units.

References:

1. Q. D. Tran, V. Sorichetti, G. Pehau-Arnaudet, M. Lenz and C. Leduc,
   "Fragmentation and entanglement limit vimentin intermediate filament
   assembly", Phys. Rev. X 13, 011014 (2023), doi:10.1103/PhysRevX.13.011014.
2. M. Panoukidou, S. Weir, V. Sorichetti, Y. Gutierrez Fosado, M. Lenz and
   D. Michieletto, "Runaway transition in irreversible polymer condensation
   with cyclization", Phys. Rev. Research 6, 023189 (2024),
   doi:10.1103/PhysRevResearch.6.023189.
"""

from .fragmentation import FragmentationParams
from .population import Population
from .rings import RingsParams

__all__ = ["FragmentationParams", "Population", "RingsParams"]
__version__ = "1.0.0"
