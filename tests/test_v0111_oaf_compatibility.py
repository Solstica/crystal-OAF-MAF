from pathlib import Path

import numpy as np

from pvdf_pf.calibration.oaf_compatibility import compatibility_series, summarize_compatibility

ROOT = Path(__file__).resolve().parents[1]
RUI = ROOT / "data/literature/rui2022"


def test_donor_progress_spans_zero_to_one():
    rows = compatibility_series(RUI, sample_state="poled")
    q = np.asarray([r.q_donor_progress for r in rows])
    assert np.isclose(q[0], 0.0, atol=1e-12)
    assert np.isclose(q[-1], 1.0, atol=1e-12)
    assert np.all(np.diff(q) >= 0.0)


def test_poled_transfer_has_irreducible_q1_endpoint_gap():
    rows = compatibility_series(RUI, sample_state="poled")
    hot = rows[-1]
    assert hot.endpoint_gap_if_q1 is not None
    # Same-state closed epsilon_OAF,eff is about 46.5 whereas Figure-S2
    # epsilon_MOAF is about 31.4 at 40 C.
    assert hot.endpoint_gap_if_q1 > 10.0
    assert hot.epsilon_roaf_required is None


def test_required_roaf_is_not_a_temperature_independent_constant():
    rows = compatibility_series(RUI, sample_state="poled")
    required = np.asarray([
        r.epsilon_roaf_required for r in rows if r.epsilon_roaf_required is not None
    ])
    assert np.all(np.isfinite(required))
    assert required.max() - required.min() > 5.0


def test_constant_low_temperature_roaf_anchor_does_not_make_transfer_compatible():
    rows = compatibility_series(RUI, sample_state="poled")
    summary = summarize_compatibility(rows)
    assert summary["status"] == "TRANSFERRED_Q_AND_FIGURE_S2_MOAF_NOT_JOINTLY_COMPATIBLE"
    assert summary["max_abs_missing_response_constant_roaf_anchor"] > 10.0
    assert summary["n_constant_anchor_convex_incompatible"] >= 1


def test_unpoled_is_also_reported_without_forcing_a_fit():
    rows = compatibility_series(RUI, sample_state="unpoled")
    summary = summarize_compatibility(rows)
    assert summary["n_points"] == 8
    assert summary["status"] == "TRANSFERRED_Q_AND_FIGURE_S2_MOAF_NOT_JOINTLY_COMPATIBLE"
