"""Direct Simulation Monte Carlo (DSMC) of Smoluchowski equations.

- :mod:`smolu_dsmc.fragmentation`: annealing and fragmentation, optionally with
  entanglement (Tran, Sorichetti et al., Phys. Rev. X 13, 011014 (2023);
  bioRxiv 10.1101/2022.03.19.484978).
- :mod:`smolu_dsmc.rings`: irreversible condensation with ring formation
  (Panoukidou et al., arXiv:2210.14010).
"""

from .fragmentation import FragmentationParams
from .population import Population
from .rings import RingsParams

__all__ = ["FragmentationParams", "Population", "RingsParams"]
__version__ = "1.0.0"
