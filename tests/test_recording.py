import json

import numpy as np
import pytest

from smolu_dsmc import FragmentationParams, fragmentation
from smolu_dsmc.cli import main
from smolu_dsmc.fragmentation import simulate
from smolu_dsmc.fragmentation_gillespie import simulate_gillespie
from smolu_dsmc.recording import TimeGrid, log_time_grid


def test_time_grid_records_in_order():
    calls = []
    grid = TimeGrid([0.0, 1.0, 2.0, 5.0], lambda i, t, pop: calls.append((i, t, pop)))
    assert grid.next_time == 0.0
    assert grid.record_until(1.5, "a") == 2.0
    assert grid.record_until(1.8, "b") == 2.0
    assert grid.record_until(4.0, "c") == 5.0
    grid.record_remaining("d")
    assert calls == [(0, 0.0, "a"), (1, 1.0, "a"), (2, 2.0, "c"), (3, 5.0, "d")]
    assert grid.next_time == float("inf")


def test_log_time_grid():
    grid = log_time_grid(1e-3, 10.0, 5)
    assert grid[0] == pytest.approx(1e-3) and grid[-1] == pytest.approx(10.0)
    assert np.allclose(np.diff(np.log10(grid)), 1.0)
    with pytest.raises(ValueError):
        log_time_grid(0.0, 10.0, 5)
    with pytest.raises(ValueError):
        log_time_grid(1.0, 10.0, 1)


@pytest.mark.parametrize("simulate_fn", [simulate, simulate_gillespie])
def test_recorded_state_matches_event_history(simulate_fn):
    params = FragmentationParams(ntot=1000, kd=1e-2, time_max=20.0, seed=5)
    events = []
    simulate_fn(params, on_event=lambda t, pop: events.append((t, params.ntot / pop.n)))
    times = log_time_grid(1e-3, 20.0, 40)
    recorded = []
    simulate_fn(
        params, record_times=times, on_record=lambda i, t, pop: recorded.append(params.ntot / pop.n)
    )
    # Same seed, same trajectory: the state at t is the one after the last event at or before t.
    event_times, lengths = np.array(events).T
    idx = np.searchsorted(event_times, times, side="right") - 1
    expected = np.where(idx >= 0, lengths[np.maximum(idx, 0)], 1.0)
    assert len(recorded) == len(times)
    assert np.array_equal(recorded, expected)


def test_cli_time_grid_and_run_record(tmp_path):
    main(
        ["frag", "--ntot", "500", "--time-max", "5", "--n-points", "12", "--outdir", str(tmp_path)]
    )
    tag = "ntot500_Kd1.00e-03_density1.00e+00"
    series = np.loadtxt(tmp_path / f"mav_vs_t_{tag}_tmax5.00e+00.dat")
    assert series.shape == (12, 3)
    assert series[0, 0] == pytest.approx(1e-5) and series[-1, 0] == pytest.approx(5.0)
    assert len(list((tmp_path / f"histograms_{tag}").iterdir())) == 12

    record = json.loads((tmp_path / f"run_{tag}_tmax5.00e+00.json").read_text())
    assert record["status"] == "finished"
    assert record["method"].startswith("dsmc")
    assert record["parameters"]["seed"] == 123 and record["parameters"]["ntot"] == 500
    assert record["recording"] == {"mode": "time_grid", "n_points": 12, "t_min": 1e-5}
    assert record["results"]["final_time"] >= 5.0
    assert record["command"][1:3] == ["frag", "--ntot"]


def test_run_record_marks_failed_runs(tmp_path, monkeypatch):
    def broken(*args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(fragmentation, "simulate", broken)
    with pytest.raises(RuntimeError):
        main(["frag", "--ntot", "100", "--time-max", "1", "--outdir", str(tmp_path)])
    (path,) = tmp_path.glob("run_*.json")
    record = json.loads(path.read_text())
    assert record["status"] == "failed" and "boom" in record["results"]["error"]


def test_rings_run_record(tmp_path):
    main(["rings", "--ntot", "100", "--outdir", str(tmp_path)])
    record = json.loads((tmp_path / "run_ntot100_density1.00e+00_pref0.50.json").read_text())
    assert record["status"] == "finished" and record["model"] == "rings"
