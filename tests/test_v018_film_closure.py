from pathlib import Path

import numpy as np

from pvdf_pf.calibration.film_closure import audit_film_closure, closure_series, summarize_closure
from pvdf_pf.calibration.source_data import load_digitized_curve_csv

ROOT = Path(__file__).resolve().parents[1]
RUI = ROOT / "data/literature/rui2022"


def test_figure3a_digitization_has_two_monotone_film_curves():
    pts = load_digitized_curve_csv(RUI / "figure_3a_film_eps_digitized.csv")
    assert len(pts) == 16
    by_sample = {}
    for sample in ("unpoled", "poled"):
        rows = sorted(
            [p for p in pts if p.sample_state == sample],
            key=lambda p: p.temperature_C,
        )
        assert len(rows) == 8
        vals = np.asarray([p.value for p in rows])
        assert np.all(np.diff(vals) > 0.0)
        by_sample[sample] = vals
    assert np.all(by_sample["poled"] > by_sample["unpoled"])


def test_project_four_component_set_still_does_not_close_measured_poled_film():
    cold = audit_film_closure(RUI, sample_state="poled", temperature_C=-30.0)
    hot = audit_film_closure(RUI, sample_state="poled", temperature_C=40.0)

    assert cold.epsilon_film_project_parallel < cold.epsilon_film_measured
    assert hot.epsilon_film_project_parallel < hot.epsilon_film_measured
    assert cold.epsilon_oaf_mean_required_for_closure > 30.0
    assert hot.epsilon_oaf_mean_required_for_closure > 60.0
    assert not cold.literal_si_available
    assert hot.literal_si_available


def test_literal_si_audit_preserves_raw_fraction_sum_and_does_not_hide_residual():
    row = audit_film_closure(RUI, sample_state="poled", temperature_C=0.0)
    assert row.literal_si_available
    assert row.literal_si_fraction_sum is not None
    assert 1.0 < row.literal_si_fraction_sum < 1.02
    assert row.epsilon_film_literal_si is not None
    assert row.residual_literal_si is not None
    assert row.residual_literal_si < -5.0


def test_project_closure_summary_flags_unresolved_gap():
    rows = closure_series(RUI, sample_state="poled")
    summary = summarize_closure(rows)
    assert summary["n_points"] == 8
    assert summary["rmse_project"] > 5.0
    assert summary["max_abs_relative_error_project"] > 0.30
    assert summary["mean_required_oaf_gap"] > 20.0
