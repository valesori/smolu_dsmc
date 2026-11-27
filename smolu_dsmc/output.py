"""Output files. Names and formats are those of the original scripts."""

from __future__ import annotations

from pathlib import Path
from typing import TextIO

import numpy as np

from .fragmentation import FragmentationParams
from .population import Population
from .rings import RingsParams, RingsState


def fragmentation_tag(params: FragmentationParams) -> str:
    tag = f"ntot{params.ntot}_Kd{params.kd:.2e}_density{params.density:.2e}"
    if params.gamma is not None:
        tag += f"_gamma{params.gamma:.2e}"
        # Runs with the legacy (reduced-density) x are marked; the default is not.
        if params.entanglement_x == "reduced":
            tag += "_legacy"
    return tag


def rings_tag(params: RingsParams) -> str:
    return f"ntot{params.ntot}_density{params.density:.2e}_pref{params.prefactor:.2f}"


class FragmentationWriter:
    """Writes the files of a fragmentation run into ``outdir``.

    ``<tag>`` is ``ntot…_Kd…_density…``, followed for entangled runs by
    ``_gamma…`` and, with the legacy convention ``entanglement_x = "reduced"``,
    ``_legacy``.

    - ``mav_vs_t_<tag>_tmax<time_max>.dat``: the waiting time, the mean mass
      and the extent of reaction ``1 - n / ntot``, either after each event
      (:meth:`on_event`) or at the times of a grid (:meth:`on_record`).
    - ``histograms_<tag>/n<index>_t<time>``: mass histograms (mass, count),
      where ``<index>`` is the Monte Carlo step (:meth:`on_snapshot`) or the
      index of the grid time (:meth:`on_record`).
    - ``masses_final_<tag>_t<time>.dat``: final masses, padded with zeros to
      ``ntot`` entries as in the original scripts.
    - ``run_<tag>_tmax<time_max>.json``: the run record, written by the CLI
      (:mod:`smolu_dsmc.provenance`).

    ``suffix`` is appended to ``<tag>``; the CLI uses ``_gillespie`` for runs of
    :mod:`smolu_dsmc.fragmentation_gillespie`, so that they are not mistaken for
    runs of the paper's DSMC algorithm.
    """

    def __init__(self, params: FragmentationParams, outdir: Path | str = ".", suffix: str = ""):
        self.params = params
        self.outdir = Path(outdir)
        tag = fragmentation_tag(params) + suffix
        self.tag = tag
        self.hist_dir = self.outdir / f"histograms_{tag}"
        self.hist_dir.mkdir(parents=True, exist_ok=True)
        self.timeseries_path = self.outdir / f"mav_vs_t_{tag}_tmax{params.time_max:.2e}.dat"
        self.run_record_path = self.outdir / f"run_{tag}_tmax{params.time_max:.2e}.json"
        self._fout: TextIO | None = None

    def __enter__(self) -> FragmentationWriter:
        self._fout = open(self.timeseries_path, "w")
        return self

    def __exit__(self, *exc) -> None:
        if self._fout is not None:
            self._fout.close()
            self._fout = None

    def on_event(self, time: float, pop: Population) -> None:
        ntot = self.params.ntot
        mean_mass = ntot / pop.n  # Mass is conserved
        extent = 1.0 - pop.n / ntot
        self._fout.write(f"{time:.4e} {mean_mass:.4e} {extent:.4e}\n")

    def on_snapshot(self, step: int, time: float, pop: Population) -> None:
        path = self.hist_dir / f"n{step:010d}_t{time:.3e}"
        np.savetxt(path, pop.histogram(), fmt="%d %d")

    def on_record(self, index: int, time: float, pop: Population) -> None:
        """Time-series line and histogram at the grid time ``time``."""
        self.on_event(time, pop)
        self.on_snapshot(index, time, pop)

    def write_final(self, time: float, pop: Population) -> Path:
        masses = np.zeros(self.params.ntot, dtype=np.int64)
        masses[: pop.n] = pop.masses()
        path = self.outdir / f"masses_final_{self.tag}_t{time:.3e}.dat"
        np.savetxt(path, masses, fmt="%d")
        return path


class RingsWriter:
    """Writes the files of a ring-formation run into ``outdir``.

    - ``mav_vs_t_<tag>.dat``: waiting time, mean mass of all molecules, of
      chains, and of rings (0 if there are none).
    - ``nchain_rings_<tag>.dat``: waiting time, number of chains and of rings.
    """

    def __init__(self, params: RingsParams, outdir: Path | str = ".", verbose: bool = False):
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        tag = rings_tag(params)
        self.mav_path = self.outdir / f"mav_vs_t_{tag}.dat"
        self.counts_path = self.outdir / f"nchain_rings_{tag}.dat"
        self.run_record_path = self.outdir / f"run_{tag}.json"
        self.verbose = verbose
        self._mav: TextIO | None = None
        self._counts: TextIO | None = None

    def __enter__(self) -> RingsWriter:
        self._mav = open(self.mav_path, "w")
        self._counts = open(self.counts_path, "w")
        return self

    def __exit__(self, *exc) -> None:
        for f in (self._mav, self._counts):
            if f is not None:
                f.close()
        self._mav = self._counts = None

    def on_event(self, time: float, state: RingsState) -> None:
        self._mav.write(
            f"{time:.4e} {state.mean_mass():.4e} "
            f"{state.mean_chain_mass():.4e} {state.mean_ring_mass():.4e}\n"
        )
        self._counts.write(f"{time:.4e} {state.n_chains} {state.n_rings}\n")
        if self.verbose:
            print(state.n_chains, state.n_rings)
