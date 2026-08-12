from pathlib import Path

import numpy as np

from pvdf_pf.calibration.project_devitrification import build_regularized_oaf_state


ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "data/literature/rui2022/figure_5b_maintext_digitized.csv"


def test_regularized_progress_endpoints():
    lo = build_regularized_oaf_state(CSV, temperature_C=-30)
    hi = build_regularized_oaf_state(CSV, temperature_C=40)
    assert np.isclose(lo.q_devitrification, 0.0)
    assert np.isclose(hi.q_devitrification, 1.0)
    assert np.isclose(lo.roaf, 0.2)
    assert np.isclose(lo.moaf, 0.0)
    assert np.isclose(hi.roaf, 0.0)
    assert np.isclose(hi.moaf, 0.2)


def test_regularized_structural_fractions_close():
    for T in (-30, -10, 0, 20, 40):
        s = build_regularized_oaf_state(CSV, temperature_C=T)
        assert np.isclose(s.crystal + s.roaf + s.moaf + s.iaf, 1.0)
        assert np.isclose(s.total_oaf, 0.2)
        assert 0.0 <= s.q_devitrification <= 1.0


def test_devitrification_progress_is_monotone_on_source_grid():
    Ts = np.arange(-30, 41, 5, dtype=float)
    q = np.asarray([build_regularized_oaf_state(CSV, temperature_C=T).q_devitrification for T in Ts])
    assert np.all(np.diff(q) >= -1e-12)


def test_raf_and_maf_progress_agree_within_digitization_uncertainty():
    for T in np.arange(-30, 41, 5, dtype=float):
        s = build_regularized_oaf_state(CSV, temperature_C=T)
        assert s.q_internal_disagreement <= 0.08


def test_outside_source_window_is_rejected():
    import pytest

    with pytest.raises(ValueError):
        build_regularized_oaf_state(CSV, temperature_C=50)
